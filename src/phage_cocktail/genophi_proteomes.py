"""Versioned, resumable gene calling for GenoPHI's locked nucleotide inputs."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import platform
import re
import time

from .artifacts import atomic_write_json, code_provenance, read_json_object, sha256_file
from .cohort import parse_dna_fasta
from .host_study import verify_resources

PYRODIGAL_VERSION = "3.7.1"
RESOURCE_LOCK = "data/studies/genophi_host_baseline.resources.v1.json"
SETS = {"host_genome": ("strain_AAs", "single"), "phage_genome": ("phage_AAs", "meta")}
CALLER_PARAMETERS = {
    "translation_table": 11,
    "closed": False,
    "mask": False,
    "min_mask": 50,
    "min_gene": 90,
    "min_edge_gene": 60,
    "max_overlap": 60,
    "backend": "detect",
    "force_nonsd": False,
    "start_weight": 4.35,
    "include_stop": False,
    "meta_bins": "Only built-in bins with training_info.translation_table == 11",
    "protein_ids": "<resource-lock id>_<1-based gene number across all contigs>",
    "fasta_width": 60,
}


def _validate_id(genome_id):
    if not isinstance(genome_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", genome_id):
        raise ValueError("genome id must be a safe FASTA/file identifier")


def call_proteome(contigs, genome_id, mode):
    """Call one genome; training never pools separate genomes or changes code 11.

    Meta mode restricts the *gene prediction* models, not just the final
    translation: Pyrodigal's default bins also include code-4 models.
    """
    import pyrodigal

    _validate_id(genome_id)
    if mode not in {"single", "meta"} or not contigs:
        raise ValueError("gene calling requires contigs and single or meta mode")
    options = {name: CALLER_PARAMETERS[name] for name in (
        "closed", "mask", "min_mask", "min_gene", "min_edge_gene", "max_overlap", "backend"
    )}
    if mode == "meta":
        bins = pyrodigal.MetagenomicBins(
            b for b in pyrodigal.METAGENOMIC_BINS if b.training_info.translation_table == 11
        )
        if not bins:
            raise ValueError("no code-11 metagenomic models available")
        finder = pyrodigal.GeneFinder(meta=True, metagenomic_bins=bins, **options)
    else:
        finder = pyrodigal.GeneFinder(meta=False, **options)
        finder.train(*contigs, translation_table=11, force_nonsd=False, start_weight=4.35)

    lines, counts, models = [], [], []
    number = 0
    for sequence in contigs:
        genes = finder.find_genes(sequence)
        if genes.training_info.translation_table != 11:
            raise ValueError("gene prediction unexpectedly used a non-11 genetic code")
        counts.append(len(genes))
        if genes.metagenomic_bin is not None:
            models.append(genes.metagenomic_bin.description)
        for gene in genes:
            protein = gene.translate(translation_table=11, include_stop=False)
            if not protein or set(protein) - set("ACDEFGHIKLMNPQRSTVWYX"):
                raise ValueError(f"invalid protein translation for {genome_id}")
            number += 1
            lines.append(f">{genome_id}_{number}\n")
            lines.extend(protein[i:i + 60] + "\n" for i in range(0, len(protein), 60))
    if not number:
        raise ValueError(f"no proteins called for {genome_id}")
    return "".join(lines).encode("ascii"), {
        "proteins": number, "contigs": len(contigs), "genes_per_contig": counts,
        "nucleotides": sum(map(len, contigs)), "meta_models_used": sorted(set(models)),
    }


def verify_proteome(path, genome_id, *, expected=None):
    """Check sequential, genome-prefixed IDs and non-empty protein sequences."""
    _validate_id(genome_id)
    count, length = 0, 0
    with Path(path).open(encoding="ascii") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count and not length:
                    raise ValueError("empty protein sequence")
                count += 1
                if line != f">{genome_id}_{count}":
                    raise ValueError("protein IDs must be genome-prefixed and sequential")
                length = 0
            elif not count or set(line) - set("ACDEFGHIKLMNPQRSTVWYX"):
                raise ValueError("invalid protein FASTA")
            else:
                length += len(line)
    if not count or not length or (expected is not None and count != expected):
        raise ValueError("empty or incomplete proteome")
    return count


def _derive_one(entry, source, output, lock_sha, caller_version):
    directory, mode = SETS[entry["kind"]]
    identity = entry["id"]
    relative = f"{directory}/{identity}.faa"
    target = output / relative
    checkpoint = output / "checkpoints" / f"{directory}.{identity}.json"
    binding = {
        "id": identity, "kind": entry["kind"], "path": relative, "mode": mode,
        "source_path": entry["path"], "source_sha256": entry["sha256"],
        "source_size_bytes": entry["size_bytes"], "resource_lock_sha256": lock_sha,
        "pyrodigal": caller_version, "parameters": CALLER_PARAMETERS,
    }
    if checkpoint.exists():
        record = read_json_object(checkpoint)
        if any(record.get(k) != v for k, v in binding.items()):
            raise ValueError(f"changed gene-calling checkpoint: {identity}")
        if (not target.is_file() or target.stat().st_size != record["size_bytes"]
                or sha256_file(target) != record["sha256"]):
            raise ValueError(f"missing or changed derived proteome: {identity}")
        verify_proteome(target, identity, expected=record["proteins"])
        return record
    if target.exists():
        raise ValueError(f"unreceipted proteome requires review: {relative}")
    started = time.perf_counter()
    contigs = parse_dna_fasta(source.read_text(encoding="ascii"), source)
    payload, stats = call_proteome(contigs, identity, mode)
    temporary = target.with_suffix(".faa.partial")
    # Exclusive creation retains stale partials for inspection after interruption.
    with temporary.open("xb") as stream:
        stream.write(payload)
    verify_proteome(temporary, identity, expected=stats["proteins"])
    temporary.rename(target)
    record = {
        **binding, **stats, "size_bytes": target.stat().st_size, "sha256": sha256_file(target),
        "wall_seconds": time.perf_counter() - started,
    }
    atomic_write_json(record, checkpoint)
    print(f"Called {identity}: {stats['proteins']} proteins ({mode}, code 11)", flush=True)
    return record


def prepare_proteomes(repo, data_root, output_dir, receipt_path, *, workers=8, source_receipt=None):
    """Verify all locked inputs, checkpoint each genome, then publish a receipt."""
    repo, output, receipt = Path(repo).resolve(), Path(output_dir).resolve(), Path(receipt_path)
    if receipt.exists():
        raise ValueError("receipt already exists; choose a new immutable receipt path")
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    caller_version = version("pyrodigal")
    if caller_version != PYRODIGAL_VERSION:
        raise ValueError(f"gene calling requires pyrodigal=={PYRODIGAL_VERSION}")
    started, utc = time.perf_counter(), datetime.now(timezone.utc).isoformat()
    lock_path = repo / RESOURCE_LOCK
    lock, lock_sha = read_json_object(lock_path), sha256_file(lock_path)
    resources = verify_resources(lock, data_root)
    entries = [r for r in lock["resources"] if r["kind"] in SETS]
    ids = [r["id"] for r in entries]
    if len(ids) != len(set(ids)):
        raise ValueError("genome IDs must be globally unique across host and phage sets")
    for identity in ids:
        _validate_id(identity)
    if source_receipt is not None:
        source_record = {
            "status": "all_locked_resources_verified", "resource_lock": RESOURCE_LOCK,
            "resource_lock_sha256": lock_sha, "verified_files": len(resources),
            "total_size_bytes": sum(r["size_bytes"] for r in lock["resources"]),
            "resources": [{k: r[k] for k in ("id", "kind", "path", "sha256", "size_bytes")}
                          for r in lock["resources"]],
        }
        if Path(source_receipt).exists():
            if read_json_object(source_receipt) != source_record:
                raise ValueError("existing source receipt differs from verified inputs")
        else:
            atomic_write_json(source_record, source_receipt)
    for directory in ["checkpoints", *[values[0] for values in SETS.values()]]:
        (output / directory).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(
            lambda entry: _derive_one(entry, resources[(entry["kind"], entry["id"])],
                                      output, lock_sha, caller_version), entries
        ))
    document = {
        "schema_version": 1, "status": "genophi_proteomes_verified",
        "resource_lock": RESOURCE_LOCK, "resource_lock_sha256": lock_sha,
        "versions": {"python": platform.python_version(), "pyrodigal": caller_version},
        "parameters": CALLER_PARAMETERS, "modes": {k: v[1] for k, v in SETS.items()},
        "started_utc": utc, "wall_seconds": time.perf_counter() - started, "workers": workers,
        "provenance": code_provenance(repo, repo / "scripts/call_genophi_genes.py"),
        "counts": {kind: {
            "genomes": sum(r["kind"] == kind for r in records),
            "proteins": sum(r["proteins"] for r in records if r["kind"] == kind),
            "size_bytes": sum(r["size_bytes"] for r in records if r["kind"] == kind),
        } for kind in SETS},
        "resources": sorted(records, key=lambda r: (r["kind"], r["id"])),
    }
    atomic_write_json(document, receipt)
    return document


def verify_proteomes(receipt, root, lock):
    """Admit exactly the locked genomes; reject extra files and changed bytes."""
    root = Path(root).resolve()
    if (receipt.get("status") != "genophi_proteomes_verified"
            or receipt.get("versions", {}).get("pyrodigal") != PYRODIGAL_VERSION
            or receipt.get("parameters") != CALLER_PARAMETERS):
        raise ValueError("unsupported proteome receipt or caller parameters")
    entries = {(r["kind"], r["id"]): r for r in lock["resources"] if r["kind"] in SETS}
    admitted = {}
    for record in receipt["resources"]:
        key = record["kind"], record["id"]
        if key in admitted or key not in entries:
            raise ValueError("duplicate or unexpected proteome identity")
        directory, mode = SETS[key[0]]
        relative = f"{directory}/{key[1]}.faa"
        source = entries[key]
        if (record["path"] != relative or record["mode"] != mode
                or record["source_sha256"] != source["sha256"]
                or record["source_size_bytes"] != source["size_bytes"]
                or record["parameters"] != CALLER_PARAMETERS
                or record["pyrodigal"] != PYRODIGAL_VERSION):
            raise ValueError("proteome source, mode or parameters differ from receipt")
        path = (root / relative).resolve()
        if (not path.is_relative_to(root) or not path.is_file()
                or path.stat().st_size != record["size_bytes"]
                or sha256_file(path) != record["sha256"]):
            raise ValueError("missing or changed derived proteome")
        verify_proteome(path, key[1], expected=record["proteins"])
        admitted[key] = path
    actual = {p.resolve() for directory, _ in SETS.values() for p in (root / directory).glob("*.faa")}
    if set(admitted) != set(entries) or actual != set(admitted.values()):
        raise ValueError("proteomes must cover exactly the locked genome sets")
    return admitted
