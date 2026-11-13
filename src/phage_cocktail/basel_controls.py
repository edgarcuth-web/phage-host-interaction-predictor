"""Sequence-bound, untrained controls for the inspected BASEL CSV reference."""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from .artifacts import canonical_json
from .basel_sequences import assembly_digest
from .cohort import parse_dna_fasta
from .host_genome import host_fourmers
from .partial_hit_bounds import bounded_hit_at_k


def vector_digest(vector):
    """Bind numerical values in the declared order, with a portable byte order."""
    return sha256(np.asarray(vector, dtype="<f8").tobytes()).hexdigest()


def feature_key(assembly, representation):
    return sha256(canonical_json(dict(assembly_sha256=assembly,
                                     representation_sha256=representation)).encode()).hexdigest()


def sequence_feature(payload, expected_assembly, representation):
    """Generate features only from the declared assembled sequence content."""
    assembly = assembly_digest(parse_dna_fasta(payload.decode("ascii")))
    if assembly != expected_assembly:
        raise ValueError("sequence differs from the declared assembly; host name cannot admit it")
    vector, qc = host_fourmers(payload)
    receipt = dict(assembly_sha256=assembly, representation_sha256=representation,
                   feature_key=feature_key(assembly, representation),
                   feature_values_sha256=vector_digest(vector), **qc)
    return vector, receipt


def validate_feature(vector, receipt, expected_assembly, representation):
    """Reject a stale assembly, different feature schema or changed vector."""
    vector = np.asarray(vector, dtype=float)
    if (receipt.get("assembly_sha256") != expected_assembly
            or receipt.get("representation_sha256") != representation
            or receipt.get("feature_key") != feature_key(expected_assembly, representation)):
        raise ValueError("feature identity differs from the requested assembly/representation")
    if (vector.ndim != 1 or not len(vector) or not np.isfinite(vector).all()
            or (vector < 0).any() or not np.isclose(vector.sum(), 1., rtol=0, atol=1e-12)
            or receipt.get("feature_values_sha256") != vector_digest(vector)):
        raise ValueError("feature vector is invalid or differs from its receipt")


def cosine_panel(host_vectors, phage_vectors):
    """Score a complete panel without labels, fitting, scaling or batched BLAS."""
    if not host_vectors or not phage_vectors:
        raise ValueError("both sequence panels must be nonempty")
    vectors = [np.asarray(v, dtype=float) for v in [*host_vectors.values(), *phage_vectors.values()]]
    size = vectors[0].size
    if any(v.ndim != 1 or len(v) != size or not np.isfinite(v).all()
           or (v < 0).any() or np.sum(v * v) <= 0 for v in vectors):
        raise ValueError("cosine requires finite nonnegative vectors of one shared dimension")
    rows = []
    for h in sorted(host_vectors):
        a = np.asarray(host_vectors[h], dtype=float)
        for p in sorted(phage_vectors):
            b = np.asarray(phage_vectors[p], dtype=float)
            score = np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b))
            rows.append(dict(bacterium=h, phage=p, score=float(score)))
    return pd.DataFrame(rows)


def expected_slot_fraction(scores, marked, k):
    """Expected fraction of top-K slots in a marked set, respecting exact ties."""
    panel = sorted(scores, key=lambda p: (-scores[p], p))
    start, selected = 0, 0.
    while start < min(k, len(panel)):
        end = start + 1
        while end < len(panel) and scores[panel[end]] == scores[panel[start]]:
            end += 1
        slots = min(k - start, end - start)
        selected += slots * sum(p in marked for p in panel[start:end]) / (end - start)
        start = end
    return selected / k


def evaluate_control(observed, scores, hosts, phages, *, ks=(1, 2, 3, 5), near_training=()):
    """Keep the full Cartesian panel and all hosts when bounding missing labels."""
    hosts, phages = list(hosts), list(phages)
    if (not hosts or not phages or len(set(hosts)) != len(hosts)
            or len(set(phages)) != len(phages)):
        raise ValueError("host and candidate rosters must be nonempty and unique")
    expected = pd.MultiIndex.from_product([hosts, phages])
    actual = pd.MultiIndex.from_frame(scores[["bacterium", "phage"]])
    known = pd.MultiIndex.from_frame(observed[["bacterium", "phage"]])
    if actual.has_duplicates or set(actual) != set(expected) or not np.isfinite(scores.score).all():
        raise ValueError("scores must cover the complete finite frozen panel")
    if known.has_duplicates or not set(known).issubset(set(expected)):
        raise ValueError("observed pairs must be unique and within the frozen panel")
    if any(pd.isna(v) or isinstance(v, (bool, np.bool_)) or v not in (0, 1) for v in observed.infection):
        raise ValueError("source labels must be observed 0/1; omit unknown pairs")
    if not ks or len(set(ks)) != len(ks) or any(type(k) is not int or not 1 <= k <= len(phages) for k in ks):
        raise ValueError("K values must be unique integers within the panel")
    if not set(near_training).issubset(phages):
        raise ValueError("near-training flags must refer to frozen candidates")
    rows = []
    for host in sorted(hosts):
        lookup = {r.phage: float(r.score) for r in scores[scores.bacterium.eq(host)].itertuples()}
        labels = {r.phage: int(r.infection) for r in observed[observed.bacterium.eq(host)].itertuples()}
        bounded = bounded_hit_at_k(phages, labels, lookup, tuple(range(1, max(ks) + 1)))
        bounds = bounded["hit_at_k_bounds"]
        for k in sorted(ks):
            rows.append(dict(
                bacterium=host, k=k, hit_lower=bounds[str(k)]["lower"], hit_upper=bounds[str(k)]["upper"],
                tests_lower=k - sum(bounds[str(j)]["upper"] for j in range(1, k)),
                tests_upper=k - sum(bounds[str(j)]["lower"] for j in range(1, k)),
                unknown_slot_fraction=expected_slot_fraction(lookup, set(phages) - set(labels), k),
                near_training_i90_slot_fraction=expected_slot_fraction(lookup, set(near_training), k),
                observed_positive=bounded["observed_positive"], observed_negative=bounded["observed_negative"],
                unknown_count=bounded["unknown_count"],
            ))
    per_host = pd.DataFrame(rows)
    metrics = ["hit_lower", "hit_upper", "tests_lower", "tests_upper", "unknown_slot_fraction",
               "near_training_i90_slot_fraction"]
    macro = {str(k): {m: float(v) for m, v in frame[metrics].mean().items()}
             for k, frame in per_host.groupby("k", sort=True)}
    diagnostics = dict(status="observed_rows_only_diagnostic_with_missingness_selection_bias",
                       average_precision=None, auroc=None)
    if set(observed.infection) == {0, 1}:
        joined = observed.merge(scores, on=["bacterium", "phage"], validate="one_to_one")
        diagnostics.update(average_precision=float(average_precision_score(joined.infection, joined.score)),
                           auroc=float(roc_auc_score(joined.infection, joined.score)))
    return dict(hosts=len(hosts), candidates=len(phages), scored_pairs=len(scores), score_coverage=1.,
                observed_pairs=len(observed), unknown_pairs=len(expected) - len(observed),
                macro=macro, observed_pair_diagnostics=diagnostics,
                uncertainty="identification_bounds_not_confidence_intervals"), per_host
