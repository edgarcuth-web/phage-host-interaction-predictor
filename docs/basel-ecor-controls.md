# BASEL–ECOR fixed controls and host assembly policy

The [v2 protocol](../data/studies/basel_ecor_evaluation_protocol.v2.json) resolves the [sequence audit](basel-ecor-sequence-audit.md) execution questions for two fixed controls. It was frozen in [commit c0c01da](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/c0c01dae8d558643ec279e8cc45387b0d8780ecf) before their implementation and first scores. Protocol v1 and its historical records remain unchanged; learned model fitting and release validation are outside this run.

| Question | Decision |
| --- | --- |
| Primary host sequence | Use the 25 publicly available WGS assemblies cited in Supplementary Table 9, bound to the audited bytes and normalized identities. |
| Original training host sequences | Run the complete 25 × 52 panel again as a prespecified assembly sensitivity. Do not select the better build after seeing scores. |
| Feature identity | Recompute the same 136 strand-canonical within-contig DNA 4-mer frequencies from each admitted sequence. Bind vectors to both assembly and representation digests; never reuse features by ECOR name. |
| Candidate roster | The pinned processed CSV defines this 52-phage reference. Table 8's supplementary-only names do not create additional observations. The reason four matrix-marked phages are absent remains unknown; this run does not reproduce a claimed 56-phage experiment. |
| Bas64 | Use verified accession MZ501081.1, 39,842 bp, agreeing with original S5. Retain Table 8's conflicting printed length without asserting its cause or changing labels. |
| Unknown outcomes | Preserve all 60 unknown pairs in the full candidate panel and bound their possible impact. |
| Near training phages at 90% | Retain all eight affected candidates and report their top-K exposure. Do not remove them to improve a metric. |

The controls are uniform random ranking (all equal scores, evaluated analytically) and host–phage 4-mer cosine similarity. The latter is a fixed composition control, not a trained interaction model or infection probability. No phenotype label enters feature generation or scoring, and no hyperparameter is selected.

Report host-macro Hit@2 identification bounds over **all 25 hosts**, with K=1,3,5 secondary, capped sequential test-cost bounds, unknown-label exposure and both assembly panels. Identification bounds describe incomplete-label uncertainty; they are not confidence intervals. Differences from random use conservative outer bounds, not a claim of statistical superiority.

Named top-2 sets use deterministic lexicographic tie breaking for inspection. Expected metrics average all within-score-tie permutations. These two outputs must not be conflated.

## Results

All three controls score the same complete 25-host × 52-phage panel, including all 60 unknown-label pairs. Each entry below is a lower–upper **identification bound**, not a confidence interval.

| Control | Host-macro Hit@2 | Expected tests, capped at 2 | Secondary Hit@5 |
| --- | ---: | ---: | ---: |
| Cosine, published WGS primary | 4.00%–4.00% | 1.960–1.960 | 24.00%–24.00% |
| Cosine, original training assembly sensitivity | 4.00%–4.00% | 1.960–1.960 | 20.00%–20.00% |
| Uniform random | 37.32%–44.06% | 1.722–1.768 | 60.10%–69.56% |

Only **1/25 hosts** has a source-positive pair in the cosine top two. No unknown label occurs in those top-two sets, so their Hit@2 bounds coincide. The conservative cosine-minus-random Hit@2 bounds are **−40.06 to −33.32 percentage points** on this inspected reference. The fixed composition control is inferior here for every compatible unknown-label assignment; it should not become the candidate recommender.

Switching host assembly changes the deterministic top-two membership for **ECOR23 only**. Both builds have the same primary Hit@2; their secondary Hit@5 differs. Retain both outcomes and the prespecified primary build.

The [report](../data/studies/results/basel_ecor_fixed_controls.v1/report.json) includes every K, cost, exposure and observed-row diagnostic. The [feature manifest](../data/studies/results/basel_ecor_fixed_controls.v1/feature_manifest.csv) and [actual 136-dimensional vectors](../data/studies/results/basel_ecor_fixed_controls.v1/sequence_features.csv) bind all **102** host/phage representations to source bytes, canonical assembly content and the frozen representation. Complete score matrices and per-host bounds are saved alongside them.

Independent [numerical QA](../data/studies/results/basel_ecor_fixed_controls.v1/NUMERIC_QA.json) verifies serialized feature hashes, recomputes **2,600 cosine scores** with scalar `math.fsum`, and checks **300 host/K rows** using a separate hypergeometric no-hit calculation. Unknowns remain absent observations. Source positivity still means published binary spot clearance; these results do not demonstrate productive infection.

## Reproduce

With project dependencies installed, no additional sequence annotation, GPU, minimap2 or PDF reader is needed for this control run:

```bash
python scripts/run_basel_controls.py \
  --fetch --data-root artifacts/basel_identity_inputs \
  --out artifacts/basel_fixed_controls_new
python scripts/verify_basel_controls.py \
  --bundle artifacts/basel_fixed_controls_new \
  --data-root artifacts/basel_identity_inputs \
  --out artifacts/basel_fixed_controls_new_QA.json
```

Omit `--fetch` for cached, byte-verified inputs. Both outputs must be new. The runner checks the frozen protocol digest directly, including in shallow checkouts or after squash merges. It exits 0 after the two allowed controls, while the older sequence-audit runner retains its historical exit-2 behavior.

The next learned comparison needs a separately frozen training procedure and actual host–phage interaction terms. Preserve these negative control results and use them as comparators; do not tune this composition heuristic on BASEL outcomes.
