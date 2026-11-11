"""Gene-calling modes, genetic codes, identities and immutable derived inputs."""
from hashlib import sha256
import json

import numpy as np
import pytest

from phage_cocktail import genophi_proteomes as gp

pyrodigal = pytest.importorskip("pyrodigal")


@pytest.fixture
def coding_dna():
    rng = np.random.default_rng(123)
    codons = ["GCT", "AAA", "GGT", "ATT", "TTC", "GAA", "CCT", "CGT"]
    return "".join("TAAAGGAGGTAAAAATG" + "".join(rng.choice(codons, 100)) + "TAA"
                   for _ in range(120))


def test_modes_train_one_assembly_and_restrict_prediction_to_code_11(monkeypatch, coding_dna):
    original = pyrodigal.GeneFinder
    configured, trained = [], []

    class RecordingFinder:
        def __init__(self, **kwargs):
            configured.append(kwargs)
            self.finder = original(**kwargs)

        def train(self, *sequences, **kwargs):
            trained.append((sequences, kwargs))
            return self.finder.train(*sequences, **kwargs)

        def find_genes(self, sequence):
            return self.finder.find_genes(sequence)

    monkeypatch.setattr(pyrodigal, "GeneFinder", RecordingFinder)
    host_contigs = [coding_dna * 2, coding_dna * 2]
    host_payload, host_stats = gp.call_proteome(host_contigs, "A", "single")
    phage_payload, phage_stats = gp.call_proteome([coding_dna], "P", "meta")
    assert len(trained) == 1
    assert trained[0][0] == tuple(host_contigs)
    assert trained[0][1]["translation_table"] == 11
    assert configured[0]["meta"] is False
    assert configured[1]["meta"] is True
    bins = configured[1]["metagenomic_bins"]
    assert bins and all(b.training_info.translation_table == 11 for b in bins)
    assert host_stats["proteins"] > 2 and phage_stats["proteins"] > 1
    assert sum(host_stats["genes_per_contig"]) == host_stats["proteins"]
    assert phage_stats["meta_models_used"]
    assert all("|11|" in name for name in phage_stats["meta_models_used"])
    host_ids = [line for line in host_payload.splitlines() if line.startswith(b">")]
    phage_ids = [line for line in phage_payload.splitlines() if line.startswith(b">")]
    assert host_ids == [f">A_{n}".encode() for n in range(1, host_stats["proteins"] + 1)]
    assert not set(host_ids) & set(phage_ids)


def test_fasta_rejects_cross_contig_number_reset_and_empty_proteins(tmp_path):
    path = tmp_path / "A.faa"
    for invalid in (">A_1\nMAA\n>A_1\nMKK\n", ">A_1\n>A_2\nMAA\n",
                    ">A_1\nMAA*\n", ">NODE_1\nMAA\n"):
        path.write_text(invalid)
        with pytest.raises(ValueError):
            gp.verify_proteome(path, "A")
    path.write_text(">A_1\nMAA\n>A_2\nMKK\n")
    assert gp.verify_proteome(path, "A", expected=2) == 2
    with pytest.raises(ValueError, match="incomplete"):
        gp.verify_proteome(path, "A", expected=3)


@pytest.fixture
def small_lock(tmp_path, coding_dna, monkeypatch):
    data = tmp_path / "raw"
    data.mkdir()
    entries = []
    for kind, identity in [("host_genome", "A"), ("phage_genome", "P")]:
        sequence = coding_dna * (3 if kind == "host_genome" else 1)
        payload = f">same_contig_name\n{sequence}\n".encode()
        (data / f"{identity}.fna").write_bytes(payload)
        entries.append({"kind": kind, "id": identity, "path": f"{identity}.fna",
                        "sha256": sha256(payload).hexdigest(), "size_bytes": len(payload)})
    lock = {"resources": entries}
    (tmp_path / "lock.json").write_text(json.dumps(lock))
    monkeypatch.setattr(gp, "RESOURCE_LOCK", "lock.json")
    monkeypatch.setattr(gp, "code_provenance", lambda *args: {})
    return data, lock


def test_receipts_bind_sources_outputs_and_resumption(tmp_path, small_lock, monkeypatch):
    data, lock = small_lock
    output, receipt = tmp_path / "derived", tmp_path / "receipt.json"
    first = gp.prepare_proteomes(tmp_path, data, output, receipt, workers=1)
    assert first["versions"]["pyrodigal"] == gp.PYRODIGAL_VERSION
    assert len(gp.verify_proteomes(first, output, lock)) == 2
    with pytest.raises(ValueError, match="already exists"):
        gp.prepare_proteomes(tmp_path, data, output, receipt, workers=1)

    def forbid_recall(*args):
        pytest.fail("verified, checkpointed genomes must not be recalled")

    monkeypatch.setattr(gp, "call_proteome", forbid_recall)
    second = gp.prepare_proteomes(tmp_path, data, output, tmp_path / "recheck.json", workers=1)
    assert first["resources"] == second["resources"]
    (output / "strain_AAs/A.faa").write_text(">A_1\nMAA\n")
    with pytest.raises(ValueError, match="changed derived"):
        gp.verify_proteomes(first, output, lock)
    with pytest.raises(ValueError, match="changed derived"):
        gp.prepare_proteomes(tmp_path, data, output, tmp_path / "third.json", workers=1)


def test_changed_raw_input_and_extra_proteomes_are_not_reused(tmp_path, small_lock):
    data, lock = small_lock
    output = tmp_path / "derived"
    record = gp.prepare_proteomes(tmp_path, data, output, tmp_path / "receipt.json", workers=1)
    (output / "phage_AAs/EXTRA.faa").write_text(">EXTRA_1\nMAA\n")
    with pytest.raises(ValueError, match="exactly"):
        gp.verify_proteomes(record, output, lock)
    (data / "A.fna").write_text(">changed\nAAAA\n")
    with pytest.raises(ValueError, match="changed locked resource"):
        gp.prepare_proteomes(tmp_path, data, output, tmp_path / "again.json", workers=1)
