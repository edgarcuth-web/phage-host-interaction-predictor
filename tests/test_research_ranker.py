from copy import deepcopy
from hashlib import sha256
from io import StringIO
from pathlib import Path
from contextlib import redirect_stderr, redirect_stdout
import json
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import canonical_content_sha256
from phage_cocktail.cohort import reverse_complement
from phage_cocktail.host_genome import CANONICAL_WORDS
from phage_cocktail.research_ranker import (
    REGISTRY, _validate_weights, candidate_rows, load_ranker, main, rank_fasta, score_feature_vector,
)

ROOT = Path(__file__).resolve().parents[1]


class ResearchRankerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "query with space.fasta"
        self.fasta.write_text(">unknown sample\nAAAANNCCCC\n>second\nTTTTACGT\n")

    def relocate(self, *, weights=True):
        target = self.root / "relocated model root"
        (target / REGISTRY).parent.mkdir(parents=True)
        shutil.copyfile(ROOT / REGISTRY, target / REGISTRY)
        if weights:
            bundle = load_ranker(ROOT, method="host-4mer")
            path = bundle["registry"]["methods"]["host-4mer"]["weights"]["path"]
            (target/path).parent.mkdir(parents=True)
            shutil.copyfile(ROOT/path,target/path)
        return target

    def test_real_frozen_model_scores_all_candidates_and_binds_output(self):
        result = rank_fasta(self.fasta,method="host-4mer")
        self.assertEqual(result["status"],"scored")
        self.assertEqual(result["coverage"]["scored_candidates"],94)
        self.assertEqual(len(result["top_candidates"]),5)
        self.assertEqual(result["content_sha256"],canonical_content_sha256(result))
        self.assertFalse(result["applicability"]["within_observed_training_length"])
        self.assertEqual(result["applicability"]["known_cohort_matches"],[])

    def test_headers_order_case_and_strand_do_not_change_sequence_features(self):
        first = rank_fasta(self.fasta,method="host-4mer")
        changed = self.root/"renamed.fasta"
        changed.write_text(">not a host identity\n"+reverse_complement("TTTTACGT").lower()+"\n>different\n"+reverse_complement("AAAANNCCCC").lower()+"\n")
        second = rank_fasta(changed,method="host-4mer")
        self.assertNotEqual(first["query"]["sequence_sha256"],second["query"]["sequence_sha256"])
        self.assertEqual(first["query"]["assembly_sha256"],second["query"]["assembly_sha256"])
        self.assertEqual(first["query"]["features"],second["query"]["features"])
        self.assertEqual(first["candidates"],second["candidates"])

    def test_scalar_fourmer_counts_do_not_cross_contigs_or_delete_ambiguities(self):
        counts={w:0 for w in CANONICAL_WORDS}
        for contig in ("AAAANNCCCC","TTTTACGT"):
            for i in range(len(contig)-3):
                word=contig[i:i+4]
                if set(word)<=set("ACGT"):
                    counts[min(word,reverse_complement(word))]+=1
        values=np.array([counts[w] for w in CANONICAL_WORDS],dtype=float)
        values/=values.sum()
        result=rank_fasta(self.fasta,method="host-4mer")
        self.assertEqual(result["query"]["qc"]["valid_windows"],7)
        self.assertEqual(result["query"]["features"]["values_sha256"],sha256(values.astype("<f8").tobytes()).hexdigest())

    def test_prevalence_is_explicitly_not_genome_conditioned(self):
        other=self.root/"other.fasta";other.write_text(">different\nACGTACGTACGTACGT\n")
        a=rank_fasta(self.fasta,method="prevalence")
        b=rank_fasta(other,method="prevalence")
        self.assertFalse(a["model"]["genome_conditioned"])
        self.assertEqual(a["candidates"],b["candidates"])

    def test_no_usable_words_is_unavailable_for_genomic_method(self):
        self.fasta.write_text(">ambiguous\nNNNRYS\n")
        result=rank_fasta(self.fasta,method="host-4mer")
        self.assertEqual(result["reason"],"no_valid_fourmer_windows")
        self.assertEqual(len(result["candidates"]),94)
        self.assertTrue(all(r["score"] is None and r["status"]=="not_estimated" for r in result["candidates"]))
        self.assertEqual(result["top_candidates"],[])
        self.assertEqual(rank_fasta(self.fasta,method="prevalence")["status"],"scored")

    def test_explicit_other_species_does_not_receive_scores(self):
        result=rank_fasta(self.fasta,method="host-4mer",host_species="Klebsiella pneumoniae")
        self.assertEqual(result["reason"],"unsupported_host_species")
        self.assertEqual(result["coverage"]["not_estimated_candidates"],94)

    def test_fasta_header_cannot_attach_an_existing_host_identity(self):
        self.fasta.write_text(">AL505\nAAAACCCC\n")
        result=rank_fasta(self.fasta,method="host-4mer")
        self.assertEqual(result["applicability"]["known_cohort_matches"],[])

    def test_relocated_bundle_reproduces_the_same_result(self):
        relocated=self.relocate()
        a=rank_fasta(self.fasta,method="host-4mer")
        b=rank_fasta(self.fasta,method="host-4mer",model_root=relocated)
        self.assertEqual(a,b)

    def test_prevalence_does_not_require_genomic_weights(self):
        relocated=self.relocate(weights=False)
        self.assertEqual(rank_fasta(self.fasta,method="prevalence",model_root=relocated)["status"],"scored")
        with self.assertRaises(OSError):load_ranker(relocated,method="host-4mer")

    def test_changed_weights_and_registry_fail_closed(self):
        relocated=self.relocate()
        bundle=load_ranker(ROOT,method="host-4mer")
        weight=relocated/bundle["registry"]["methods"]["host-4mer"]["weights"]["path"]
        weight.write_bytes(weight.read_bytes()+b"\n")
        with self.assertRaisesRegex(ValueError,"pinned bytes"):load_ranker(relocated,method="host-4mer")
        registry=relocated/REGISTRY;registry.write_bytes(registry.read_bytes()+b"\n")
        with self.assertRaisesRegex(ValueError,"registry differs"):load_ranker(relocated,method="prevalence")

    def test_unknown_candidate_and_feature_order_are_rejected(self):
        bundle=load_ranker(ROOT,method="host-4mer")
        model=deepcopy(bundle["weights"])
        key=next(iter(model["phages"]));model["phages"]["unknown_phage"]=model["phages"].pop(key)
        with self.assertRaisesRegex(ValueError,"exact frozen candidate"):_validate_weights(model,bundle["registry"])
        model=deepcopy(bundle["weights"]);model["preprocessing"]["numeric_columns"].reverse()
        with self.assertRaisesRegex(ValueError,"feature schema/order"):_validate_weights(model,bundle["registry"])

    def test_partial_nonfinite_scores_and_missing_features_are_rejected(self):
        bundle=load_ranker(ROOT,method="prevalence")
        scores=score_feature_vector(bundle);scores.pop(next(iter(scores)))
        with self.assertRaisesRegex(ValueError,"complete finite"):candidate_rows(bundle["registry"]["library"],scores,5)
        scores=score_feature_vector(bundle);scores[next(iter(scores))]=float("nan")
        with self.assertRaisesRegex(ValueError,"complete finite"):candidate_rows(bundle["registry"]["library"],scores,5)
        with self.assertRaisesRegex(ValueError,"136-value"):score_feature_vector(load_ranker(ROOT,method="host-4mer"),np.ones(135)/135)

    def test_exact_ties_retain_rank_bounds_and_expected_slots(self):
        bundle=load_ranker(ROOT,method="prevalence")
        scores={r["phage_id"]:0. for r in bundle["registry"]["library"]["candidates"]}
        rows=candidate_rows(bundle["registry"]["library"],scores,5)
        self.assertTrue(all(r["rank_min"]==1 and r["rank_max"]==94 for r in rows))
        self.assertAlmostEqual(sum(r["expected_top_k_slot"] for r in rows),5.)
        self.assertEqual([r["phage_id"] for r in rows],sorted(scores))

    def test_saved_al505_vector_reproduces_all_frozen_scores(self):
        path=ROOT/'data/studies/results/paired_interaction.v1/training_features.csv.gz'
        vectors=pd.read_csv(path,float_precision="round_trip").set_index(['entity_kind','entity_id'])
        vector=vectors.loc[('host_genome','AL505'),CANONICAL_WORDS].to_numpy(float)
        saved=pd.read_csv(ROOT/'data/studies/results/new_host_retrieval.v1/scores.csv.gz',float_precision="round_trip")
        for method,source in [('host-4mer','host_4mer'),('prevalence','phage_prevalence')]:
            expected=saved[saved.panel.eq('grouped_host')&saved.method.eq(source)&saved.bacterium.eq('AL505')].set_index('phage').score
            scores=score_feature_vector(load_ranker(ROOT,method=method),vector)
            self.assertEqual(len(expected),94)
            self.assertLessEqual(max(abs(scores[p]-v) for p,v in expected.items()),1e-12)

    def test_cli_writes_new_json_and_preserves_existing_output(self):
        out=self.root/'ranking.json'
        args=['--fasta',str(self.fasta),'--method','host-4mer','--out',str(out)]
        with redirect_stdout(StringIO()):self.assertEqual(main(args),0)
        original=out.read_bytes()
        with redirect_stderr(StringIO()):self.assertEqual(main(args),2)
        self.assertEqual(out.read_bytes(),original)

    def test_invalid_dna_fastq_and_top_k_are_rejected(self):
        for contents in ('@read\nACGT\n+\nIIII\n','>protein\nMKWVTFISLLFLFSSAYS\n','>empty\n'):
            self.fasta.write_text(contents)
            with self.assertRaises(ValueError):rank_fasta(self.fasta,method="host-4mer")
        self.fasta.write_text('>valid\nACGTACGT\n')
        with self.assertRaisesRegex(ValueError,"top-k"):rank_fasta(self.fasta,method="host-4mer",top_k=95)


if __name__ == "__main__":unittest.main()
