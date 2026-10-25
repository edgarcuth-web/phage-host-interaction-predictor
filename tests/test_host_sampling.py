"""Nested sampling support, count weighting and fixed-panel integration."""
import json

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.artifacts import sha256_file
from phage_cocktail.cohort import reverse_complement
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.host_sampling import sampled_fourmers
from phage_cocktail.retrieval_study import run_retrieval_study
from test_host_genome import genome_study
from test_retrieval_study import local_study


def test_nested_sampling_counts_valid_windows_without_inventing_adjacency():
    payload = b">one\nAAAANNNNAAAA\n>two\nCCCCGGGGCCCC\n>short\nAT\n"
    vectors, audit, windows = sampled_fourmers(payload, window_bp=4, pool_windows=32, budgets=[1,4,32])
    a, c = CANONICAL_WORDS.index("AAAA"), CANONICAL_WORDS.index("CCCC")
    assert vectors[1][a] == 1
    np.testing.assert_allclose(vectors[4][[a,c]], [1/3,2/3])
    np.testing.assert_allclose(vectors[32][[a,c]], [2/5,3/5])
    assert audit.set_index("window_budget").loc[32,"selected_windows"] == 6
    assert audit.set_index("window_budget").loc[4,"valid_fourmer_windows"] == 3
    assert audit.set_index("window_budget").loc[4,"ambiguous_bases"] == 4
    assert audit.assembly_bases.eq(26).all()
    sets = {k:set(v.pool_ordinal) for k,v in windows.groupby("window_budget")}
    assert sets[1] <= sets[4] <= sets[32]
    assert windows.length.eq(4).all()


def test_sampling_is_independent_of_contig_names_order_and_strand():
    contigs = ["AACCGGTTAACCNNNTT", "CGTTTTAACCGG"]
    a = (">a\n" + contigs[0] + "\n>b\n" + contigs[1] + "\n").encode()
    b = (">x\n" + reverse_complement(contigs[1]) + "\n>y\n" + reverse_complement(contigs[0]) + "\n").encode()
    left = sampled_fourmers(a,window_bp=4,pool_windows=32,budgets=[1,4,32])
    right = sampled_fourmers(b,window_bp=4,pool_windows=32,budgets=[1,4,32])
    for k in left[0]:
        np.testing.assert_array_equal(left[0][k],right[0][k])
    pd.testing.assert_frame_equal(left[1],right[1])
    pd.testing.assert_frame_equal(left[2],right[2])
    with pytest.raises(ValueError, match="no full-length"):
        sampled_fourmers(b">a\nACG\n",window_bp=4,pool_windows=32,budgets=[1])
    with pytest.raises(ValueError, match="no valid"):
        sampled_fourmers(b">a\nNNNN\n",window_bp=4,pool_windows=32,budgets=[1])


def test_sampling_runner_reproduces_old_models_and_keeps_all_panels(genome_study):
    repo,source,output,plan=genome_study
    old=run_retrieval_study(repo,source,output,genome=True)
    parent=repo/'data/studies/genome_host_retrieval.v1.json'
    amendment={'parent_protocol':str(parent.relative_to(repo)), 'parent_sha256':sha256_file(parent),
               'sampling_control':{'window_bp':4,'pool_windows':32,'budgets':[1,4,32]},
               'overrides':{'study_id':'fixture-sampling'}}
    (repo/'data/studies/host_window_sampling.v1.json').write_text(json.dumps(amendment))
    new=output.with_name('sampling')
    report=run_retrieval_study(repo,source,new,sampling=True)
    historical=[r for r in report['metrics'] if not r['model'].startswith('host_window')]
    assert historical==old['metrics']
    assert len(report['metrics'])==63
    assert report['audit']['sampling_control_rows']==27
    scored=pd.read_csv(new/'predictions.csv.gz')
    assert set(scored.groupby(['scenario','partition','model','bacterium']).size())=={3}
    for name,digest in report['output_sha256'].items():
        assert sha256_file(new/name)==digest
    pd.testing.assert_frame_equal(pd.read_csv(output/'host_splits.csv'),pd.read_csv(new/'host_splits.csv'))


def test_sampling_parent_tampering_is_rejected_before_outputs(genome_study):
    repo,source,output,plan=genome_study
    parent=repo/'data/studies/genome_host_retrieval.v1.json'
    amendment={'parent_protocol':str(parent.relative_to(repo)), 'parent_sha256':sha256_file(parent),
               'sampling_control':{'window_bp':4,'pool_windows':32,'budgets':[1,4,32]},'overrides':{}}
    (repo/'data/studies/host_window_sampling.v1.json').write_text(json.dumps(amendment))
    parent.write_text(parent.read_text()+'\n')
    with pytest.raises(ValueError,match='parent protocol'):
        run_retrieval_study(repo,source,output,sampling=True)
    assert not output.exists()
