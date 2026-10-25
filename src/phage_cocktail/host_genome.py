"""Fixed, label-free host genome composition for known-phage retrieval."""
from itertools import product

import numpy as np
import pandas as pd

from .cohort import canonical_assembly_digest, parse_dna_fasta, reverse_complement

WORDS = list(map("".join, product("ACGT", repeat=4)))
CANONICAL_WORDS = sorted({min(w, reverse_complement(w)) for w in WORDS})
WORD_BINS = np.array([CANONICAL_WORDS.index(min(w, reverse_complement(w))) for w in WORDS])
BASE_CODES = np.full(256, 4, dtype=np.uint16)
BASE_CODES[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)


def host_fourmers(payload):
    """Exact strand-canonical frequencies, without joining contigs or deleting Ns.

    Vectorized counting agrees with the existing scalar phage implementation.
    Ambiguous windows are skipped; short contigs contribute no windows.
    """
    contigs = parse_dna_fasta(payload.decode("ascii"))
    counts = np.zeros(len(CANONICAL_WORDS), dtype=np.int64)
    bases = windows = ambiguous = 0
    for contig in contigs:
        encoded = BASE_CODES[np.frombuffer(contig.encode("ascii"), dtype=np.uint8)]
        bases += len(encoded)
        ambiguous += int((encoded == 4).sum())
        windows += max(0, len(encoded) - 3)
        if len(encoded) < 4:
            continue
        valid = ((encoded[:-3] < 4) & (encoded[1:-2] < 4)
                 & (encoded[2:-1] < 4) & (encoded[3:] < 4))
        indices = (64 * encoded[:-3] + 16 * encoded[1:-2]
                   + 4 * encoded[2:-1] + encoded[3:])[valid]
        counts += np.bincount(WORD_BINS[indices], minlength=len(counts))
    total = int(counts.sum())
    if not total:
        raise ValueError("host genome has no valid 4-mer windows")
    return counts / total, {"contigs": len(contigs), "bases": bases,
                            "ambiguous_bases": ambiguous, "possible_windows": windows,
                            "valid_windows": total}


def host_genome_features(manifest, resources):
    """Bind each genome to the previously frozen assembly identity before encoding."""
    vectors, audit = [], []
    hosts = manifest.sort_values("bacterium")
    for row in hosts.itertuples(index=False):
        path = resources[("host_genome", row.bacterium)]
        digest = canonical_assembly_digest(path)
        if digest != row.assembly_sha256:
            raise ValueError(f"host assembly differs from frozen manifest: {row.bacterium}")
        vector, qc = host_fourmers(path.read_bytes())
        vectors.append(vector)
        audit.append({"bacterium": row.bacterium, "assembly_sha256": digest, **qc})
    return (pd.DataFrame(vectors, index=hosts.bacterium,
                         columns=["host_4mer:" + word for word in CANONICAL_WORDS]),
            pd.DataFrame(audit))
