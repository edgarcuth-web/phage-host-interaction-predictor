# BASEL–ECOR fixed controls and host assembly policy

The [v2 protocol](../data/studies/basel_ecor_evaluation_protocol.v2.json) resolves the [sequence audit](basel-ecor-sequence-audit.md) execution questions for two fixed controls. It is frozen before their first scores. Protocol v1 and its historical records remain unchanged; learned model fitting and release validation are outside this run.

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
