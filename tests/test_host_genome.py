"""Genome feature semantics, frozen assembly identity and retrieval integration."""
import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import sha256_file
from phage_cocktail.cohort import canonical_assembly_digest, reverse_complement
from phage_cocktail.host_genome import host_fourmers, host_genome_features
from phage_cocktail.paired_baseline import canonical_fourmers
from phage_cocktail.retrieval_study import run_retrieval_study, score_known_panel
from test_retrieval_study import local_study  # reuse the complete table-only fixture

ROOT = Path(__file__).resolve().parents[1]


def test_vectorized_counter_matches_independent_scalar_implementation():
    # Includes all possible 4-mers and random ambiguity codes/contig boundaries.
    rng = np.random.default_rng(42)
    payloads = [(f">w\n{''.join(w)}\n").encode() for w in product("ACGT", repeat=4)]
    payloads += [(f">a\n{''.join(rng.choice(list('ACGTRYSWKMBDHVN'), 1000))}\n"
                  ">b\nACGTACGT\n>short\nT\n").encode() for _ in range(10)]
    for payload in payloads:
        vector, _ = host_fourmers(payload)
        np.testing.assert_array_equal(vector, canonical_fourmers(payload))


def test_ambiguities_are_not_deleted_and_contigs_are_not_joined():
    vector, qc = host_fourmers(b">a\naaaNttt\n>b\nacgt\n>c\na\n")
    reference, _ = host_fourmers(b">only\nACGT\n")
    np.testing.assert_array_equal(vector, reference)
    assert qc == dict(contigs=3, bases=12, ambiguous_bases=1,
                      possible_windows=5, valid_windows=1)
    with pytest.raises(ValueError, match="no valid"):
        host_fourmers(b">a\nAAA\n>b\nTTT\n")
    with pytest.raises(ValueError, match="invalid"):
        host_fourmers(b">a\nACGT-ACGT\n")


def test_counter_is_invariant_to_strand_and_contig_order():
    contigs = ["ACGTNAACCGGT", "TATTAACGTG", "AA"]
    left = ('>a\n' + '\n>b\n'.join(contigs)).encode()
    right = ('>z\n' + '\n>x\n'.join(reverse_complement(c) for c in reversed(contigs))).encode()
    a, aq = host_fourmers(left)
    b, bq = host_fourmers(right)
    np.testing.assert_array_equal(a, b)
    assert aq == bq


def test_genome_features_check_frozen_assembly_and_stable_order(tmp_path):
    paths = {}
    rows = []
    for name in ['B', 'A']:
        path = tmp_path / f'{name}.fna'
        path.write_bytes(b'>test\nAACCGGTTAACCGGTT\n')
        paths['host_genome', name] = path
        rows.append(dict(bacterium=name, assembly_sha256=canonical_assembly_digest(path)))
    manifest = pd.DataFrame(rows)
    features, audit = host_genome_features(manifest, paths)
    assert features.index.tolist() == ['A', 'B']
    assert features.shape == (2, 136)
    np.testing.assert_array_equal(features.loc['A'], features.loc['B'])
    assert audit.bacterium.tolist() == ['A', 'B']
    paths['host_genome', 'A'].write_bytes(b'>changed\nTTTTTTTT\n')
    with pytest.raises(ValueError, match='differs from frozen'):
        host_genome_features(manifest, paths)


def test_sequence_standardization_uses_only_training_hosts():
    train = pd.DataFrame(dict(bacterium=list('ABCD'), phage=['P']*4, infection=[0,0,1,1]))
    candidates = pd.DataFrame(dict(bacterium=['E'], phage=['P']))
    f = pd.DataFrame({'host_4mer:A':[0., .1, .2, .3, 1.]}, index=list('ABCDE'))
    m, s = score_known_panel(train, candidates, f, list(f), [], dict(C=.1, max_iter=2000))
    changed = f.copy()
    changed.loc['E'] = 100
    other, _ = score_known_panel(train, candidates, changed, list(f), [], dict(C=.1, max_iter=2000))
    assert m == other
    assert m['preprocessing']['means'] == [pytest.approx(.15)]
    _, labeled = score_known_panel(train, candidates.assign(infection=1), f, list(f), [],
                                   dict(C=.1, max_iter=2000))
    np.testing.assert_array_equal(s['host_logistic'], labeled['host_logistic'])


@pytest.fixture
def genome_study(local_study):
    repo, source, output, plan = local_study
    manifest_path = repo / plan['host_manifest']
    manifest = pd.read_csv(manifest_path)
    digests = []
    for i, name in enumerate(manifest.bacterium):
        path = source / f'{name}.fna'
        path.write_text(f'>host\n{"ACGT" * (i+1)}{"AAAA" * (10-i)}\n')
        digests.append(canonical_assembly_digest(path))
        plan['resources'].append(dict(kind='host_genome', id=name, path=path.name,
                                      sha256=sha256_file(path), size_bytes=path.stat().st_size))
    manifest['assembly_sha256'] = digests
    manifest.to_csv(manifest_path, index=False)
    plan['repo_files_sha256'][plan['host_manifest']] = sha256_file(manifest_path)
    (repo / 'data/studies/genome_host_retrieval.v1.json').write_text(json.dumps(plan))
    return repo, source, output, plan


def test_genome_runner_keeps_identical_panels_and_denominators(genome_study):
    repo, source, output, _ = genome_study
    report = run_retrieval_study(repo, source, output, genome=True)
    assert report['audit']['host_genomes_verified'] == 9
    assert report['audit']['genome_features'] == 136
    assert len(report['metrics']) == 42  # seven panels, six methods
    predictions = pd.read_csv(output / 'predictions.csv.gz')
    assert set(predictions.groupby(['scenario','partition','model','bacterium']).size()) == {3}
    assert len(list(output.glob('model_*.json'))) == 15
    for name, digest in report['output_sha256'].items():
        assert sha256_file(output / name) == digest
    assert len(report['contrasts']) == 7 * 6 * 2
    for scenario in {m['scenario'] for m in report['metrics']}:
        for partition in ('validation','test'):
            rows = [m for m in report['metrics'] if m['scenario']==scenario and m['partition']==partition]
            if rows:
                assert len({m['complete_panel_hosts'] for m in rows}) == 1


def test_tampered_genome_cannot_publish_a_run(genome_study):
    repo, source, output, _ = genome_study
    (source / 'A.fna').write_bytes(b'>changed\nAAAA\n')
    with pytest.raises(ValueError, match='missing or changed'):
        run_retrieval_study(repo, source, output, genome=True)
    assert not output.exists()
