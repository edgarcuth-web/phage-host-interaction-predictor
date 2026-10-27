"""Label-free, nested-window 4-mer controls for representation sampling cost."""
from hashlib import sha256

import numpy as np
import pandas as pd

from .cohort import canonical_assembly_digest, parse_dna_fasta
from .evo2_host import host_windows
from .host_genome import CANONICAL_WORDS, host_fourmers


def sampled_fourmers(payload, *, window_bp, pool_windows, budgets):
    """Select full windows from the fixed pool and count without crossing edges.

    Budget 1 uses the first full window, matching the hosted engineering pilot.
    Larger budgets are evenly spaced within the same full-window pool. Budgets
    1,4,32 are nested when the pool has at most 32 members. Short sampled tails
    are not padded or concatenated; actual support is reported for every budget.
    """
    if window_bp < 4 or not budgets or any(type(k) is not int or not 1 <= k <= pool_windows for k in budgets):
        raise ValueError("sampling requires window length >=4 and valid positive budgets")
    if len(set(budgets)) != len(budgets):
        raise ValueError("sampling budgets must be unique")
    windows = host_windows(payload, window_bp=window_bp, max_windows_per_host=pool_windows)
    full = [(i, c, start, s) for i, (c, start, s) in enumerate(windows) if len(s) == window_bp]
    if not full:
        raise ValueError("host has no full-length sampled window")
    assembly_bases = sum(map(len, parse_dna_fasta(payload.decode("ascii"))))
    vectors, audit, coordinates = {}, [], []
    for budget in budgets:
        count = min(budget, len(full))
        indices = [0] if count == 1 else [i * (len(full) - 1) // (count - 1) for i in range(count)]
        selected = [full[i] for i in indices]
        fasta = "".join(f">window_{n}\n{s}\n" for n, (_, _, _, s) in enumerate(selected)).encode("ascii")
        # Each sampled window is its own contig. Normalization weights valid
        # 4-mer counts, not equal per-window means or artificial concatenations.
        vector, qc = host_fourmers(fasta)
        vectors[budget] = vector
        audit.append({"window_budget": budget, "sampled_pool_windows": len(windows),
                      "full_pool_windows": len(full), "selected_windows": count,
                      "selected_bases": qc["bases"], "assembly_bases": assembly_bases,
                      "sampled_assembly_fraction": qc["bases"] / assembly_bases,
                      "ambiguous_bases": qc["ambiguous_bases"], "valid_fourmer_windows": qc["valid_windows"]})
        coordinates.extend({"window_budget": budget, "pool_ordinal": ordinal, "contig_ordinal": contig,
                            "start": start, "length": len(s), "input_sha256": sha256(s.encode()).hexdigest()}
                           for ordinal, contig, start, s in selected)
    return vectors, pd.DataFrame(audit), pd.DataFrame(coordinates)


def sampling_features(manifest, resources, settings):
    hosts = manifest.sort_values("bacterium")
    matrices = {k: [] for k in settings["budgets"]}
    audit, coordinates = [], []
    for row in hosts.itertuples(index=False):
        path = resources[("host_genome", row.bacterium)]
        if canonical_assembly_digest(path) != row.assembly_sha256:
            raise ValueError("sampling host assembly differs from frozen manifest")
        vectors, qc, windows = sampled_fourmers(path.read_bytes(), **settings)
        for budget, vector in vectors.items():
            matrices[budget].append(vector)
        audit.append(qc.assign(bacterium=row.bacterium, assembly_sha256=row.assembly_sha256))
        coordinates.append(windows.assign(bacterium=row.bacterium))
    features = {k: pd.DataFrame(np.stack(values), index=hosts.bacterium,
                                columns=["host_4mer:" + word for word in CANONICAL_WORDS])
                for k, values in matrices.items()}
    return features, pd.concat(audit, ignore_index=True), pd.concat(coordinates, ignore_index=True)
