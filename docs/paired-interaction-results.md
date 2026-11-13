# Explicit host–phage interaction comparison

**The primary improvement criterion was not met.** The bilinear model improves BASEL candidate ranking over the new additive model, but has lower cold-phage development AP and a gain interval crossing zero. Retain both findings; this is a research model, not a promoted recommender.

## Model and fixed training procedure

The [protocol](../data/studies/paired_interaction.v1.json) was frozen in [commit 1b5c7cb](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/1b5c7cb314845067ef944231259cd3c36e659e4b) before implementation and first fit. All three models use the same 136 host and 136 phage within-contig, strand-canonical DNA 4-mer frequencies:

| Model | Raw logit | Parameters |
| --- | --- | ---: |
| Phage only | `c + bᵀp` | 137 |
| Additive | `c + aᵀh + bᵀp` | 273 |
| Bilinear | `c + aᵀh + bᵀp + hᵀWp` | 18,769 |

The bilinear interaction changes the phage ranking with the host. It is a learned composition model; its terms do not identify receptors, defense systems or causal infection mechanisms. Training uses exact matrix-free gradients rather than storing a dense observation-by-18,496 outer-product table.

- **Cohort:** 177 original hosts, 94 phage IDs / 89 unique phage assemblies, 16,597 reconciled observed MLC labels. All 41 unknown MLC readings are excluded with reasons. Duplicate-genome label disagreements remain separate source observations.
- **Split:** reuse the frozen 95%-identity / seed-42 crossed manifest. Fit 6,834 rows; select only on the 1,746 cold-phage validation rows. No seed search.
- **Features:** fit separate host/phage means and population scales on unique fit assemblies: 117 hosts and 61 phage assemblies for development. No validation, test or BASEL scaling statistics.
- **Selection:** L2 grid 0.1, 0.01, 0.001; choose validation weighted raw-logit AP, with stronger regularization on ties. Selected L2 is 0.001 for phage-only/additive and 0.01 for bilinear under both weightings.
- **Weights:** original observations equally weighted primary; inverse exact-sequence-pair multiplicity sensitivity. Both use the same protocol and are reported.
- **Full refit:** retain selected L2 values, refit all three models on all original observed rows, and recompute full-cohort scalers from 177 host / 89 phage assemblies. Generate every BASEL score matrix before its labels enter evaluation.

All 24 fits converged under the frozen stopping/gradient criteria. Twelve selected development/full-refit models are saved. Raw logits determine ranks; sigmoid outputs are uncalibrated diagnostics only.

## Development test results

Primary observation-equal AP, using identical held-out observations for each model:

| Regime | Phage only | Additive | Bilinear |
| --- | ---: | ---: | ---: |
| **Cold phage, primary** | 0.4630 | **0.5516** | 0.5063 |
| Cold host | 0.4786 | **0.5005** | 0.4586 |
| Double cold | 0.3692 | **0.4479** | 0.3776 |
| Warm pair | 0.4525 | **0.7157** | 0.5754 |

The primary bilinear-minus-additive AP difference is **−0.0454**, with conditional 95% phage-group bootstrap interval **[−0.1084, 0.0979]** over ten test phage groups. The sequence-pair-equal sensitivity gives **−0.0520 [−0.1124, 0.0956]**. Neither lower bound exceeds zero, so improvement is **not established**. These are conditional development intervals on already inspected data, not independent release evidence.

## Published BASEL reference

All 25 hosts × 52 phages remain scored, with 60 source-unknown pairs retained. Primary host assemblies are the published WGS sets from the [v2 assembly policy](basel-ecor-controls.md).

| Primary observation-equal model/control | Macro Hit@2 identification bounds | Expected tests capped at two |
| --- | ---: | ---: |
| Phage-only learned | 12%–12% | 1.880–1.880 |
| Additive learned | 12%–12% | 1.880–1.880 |
| **Bilinear learned** | **48%–60%** | **1.76–1.80** |
| Frozen uniform random | 37.32%–44.06% | 1.722–1.768 |
| Frozen cosine | 4%–4% | 1.960–1.960 |

These are missing-label identification bounds, **not confidence intervals**. The primary bilinear-minus-additive bounds are **+36 to +48 percentage points**, but this favorable reference outcome does not override the failed development criterion.

The bilinear training-assembly sensitivity is **48%–68%**. Sequence-pair-equal training gives **48%–64%** under both host builds. The public-WGS build remains primary. No assembly, L2, seed or method is selected from these outcomes.

For the primary bilinear ranking, 10% of top-two slots have unknown labels and 28% belong to the eight phages connected to training groups at the 90% identity sensitivity. Retain all 52; no favorable subset is substituted. Training targets MLC positivity, whereas BASEL labels describe binary spot clearance. The endpoints and sample composition differ, so the two results are not interchangeable or proof of productive infection.

## Weights, verification and reproduction

The [result bundle](../data/studies/results/paired_interaction.v1/report.json) contains actual training vectors and feature identities, source observations/exclusions, all validation trials, metrics, full score matrices and both weightings/builds. [models.json](../data/studies/results/paired_interaction.v1/models.json) names each encoder and coefficient array in [weights.npz](../data/studies/results/paired_interaction.v1/weights.npz). Load NumPy arrays with `allow_pickle=False`; no Python objects are serialized.

Independent [numerical QA](../data/studies/results/paired_interaction.v1/NUMERIC_QA.json) verifies 271 feature vectors, all twelve models' observation-level gradients/objectives, 53,334 development and 15,600 BASEL scores, validation-only parameter selection, primary bootstrap decisions and 1,200 BASEL host/K rows. Matrix products and sklearn AP/AUROC independently check the production scorer and AP implementation. The existing historical studies, weights and protocols are unchanged.

Use **Python 3.12.14** and the exact numerical package versions in the frozen protocol. The real-data study is CPU-only with one BLAS thread:

```bash
python scripts/run_paired_interaction.py \
  --fetch --data-root artifacts/interaction_inputs \
  --out artifacts/paired_interaction_new
python scripts/verify_paired_interaction.py \
  --bundle artifacts/paired_interaction_new \
  --data-root artifacts/interaction_inputs \
  --out artifacts/paired_interaction_new_QA.json
```

Omit `--fetch` for cached verified inputs; outputs must be new. Missing sequences, changed feature schemas, nonconvergence or altered frozen inputs fail instead of producing a reduced comparison. No serving adapter or SDL admission is added.

The next experiment should change a scientifically justified feature or data limitation under a new protocol. Keep this baseline and its unfavorable development result; do not tune against the inspected BASEL labels or treat the positive reference ranking alone as validation.
