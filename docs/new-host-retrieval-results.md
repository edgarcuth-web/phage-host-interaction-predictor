# Matched new-host, 94-phage comparison

**The bilinear model did not establish the frozen hit or cost improvement.**
It performs below train-only phage prevalence on both test panels. Keep the
negative result and the stronger existing comparators; no model is promoted.

The [protocol](../data/studies/new_host_retrieval.v1.json) was frozen in PR #47.
The [implementation](../data/studies/new_host_retrieval.implementation.v1.json)
was committed as [6377eda](https://github.com/edgarcuth-web/phage-host-interaction-predictor/commit/6377eda61c50350a66995f6371e97080b3f026c3)
before the first real-data fit or score. Historical studies remain separate.

## Test results

Every method scores the same 94 named candidates. Primary metrics include all
34 complete grouped-test hosts and all 61 complete B2-test hosts, including one
all-negative host in each panel. B2's ECOR53 is scored but excluded whole from
complete-panel metrics because two labels are unknown.

| Method | Primary Hit@2 | Primary tests capped at 2 | B2 Hit@2 | B2 tests capped at 2 |
| --- | ---: | ---: | ---: | ---: |
| Uniform random expectation | 43.71% | 1.734 | 48.31% | 1.700 |
| Train-only phage prevalence | 73.53% | 1.441 | **90.16%** | 1.254 |
| Annotation proxy | 82.35% | 1.353 | 80.33% | 1.344 |
| Host 4-mer | 88.24% | 1.235 | 75.41% | 1.410 |
| GenoPHI, frozen scores | **91.18%** | **1.176** | 78.69% | 1.295 |
| Phage-only sequence model | 55.88% | 1.500 | 85.25% | **1.230** |
| Additive sequence model | 55.88% | 1.500 | 85.25% | **1.230** |
| Bilinear sequence model | 64.71% | 1.412 | 63.93% | 1.443 |

One test is one equal-cost unit. Under the capped stop-at-first-positive policy,
unsuccessful hosts still consume the full budget. Hit probability and consumption
are distinct outcomes; the B2 phage-only control illustrates that tradeoff.

| Bilinear minus prevalence | Primary grouped test | B2 stress test |
| --- | --- | --- |
| Hit@2 difference, percentage points | −8.82 [−25.00, +6.90] | −26.23 [−38.89, −13.46] |
| Two-test consumption difference | −0.029 [−0.143, +0.067] | +0.189 [+0.076, +0.290] |

Brackets are conditional 95% paired host-group intervals with 2,000 draws.
They describe inspected development data, not independent release validation.
Neither primary gate passes; the joint hit-and-cost benefit is not established.
Validation is also retained: bilinear Hit@2 is 73.08% versus prevalence's 84.62%.
All K, eligible-only summaries, intervals and subgroup failures are in the
[result bundle](../data/studies/results/new_host_retrieval.v1/).

## Training and failure analysis

- All 21 pooled fits converged. Primary L2 is 0.1 for all three kinds; B2 selects
  0.001 for phage-only/additive and 0.1 for bilinear using only its frozen inner
  94-fit/21-validation-host split, then refits on 115 outer-training hosts.
- The 282 per-phage estimator slots use 277 logistic fits and five single-class
  prevalence fallbacks. Primary annotation weights and all GenoPHI scores are
  reused; 271 sequence vectors pass their frozen identities.
- Bilinear improves Hit@2 over prevalence for two primary-test hosts, loses for
  five and ties for 27. On B2 it loses for 16, improves for none and ties for 45.
  The primary net Hit@2 deficit is in the B2 and D subgroups; small counts do not
  establish a biological cause.
- A concrete follow-up hypothesis is to preserve the measured train-only
  per-phage prior and test a host-dependent residual against it. The phage-only
  sequence control already loses ranking information versus prevalence on the
  primary panel. This motivates a separately frozen experiment; it does not
  prove that a residual will help or identify the physical source of failure.

No additional seed, feature set or test-selected L2 is substituted. The next
planned deliverable is the research FASTA-to-ranking command with declared
baselines and applicability, alongside public new-host data intake. The bilinear
model is not selected as its default recommender.

## Artifacts and reproduction

[models.json](../data/studies/results/new_host_retrieval.v1/models.json) binds all
21 pooled coefficient arrays, the three fit-only encoders, selected models and
the three newly fitted per-phage JSON models. Load
[weights.npz](../data/studies/results/new_host_retrieval.v1/weights.npz) with
`allow_pickle=False`. Complete target and selection scores, every trial,
per-host/candidate diagnostics, support and exclusions are retained. Named
top-five lists use lexicographic ties for inspection; metrics average exact ties.

Independent [numerical QA](../data/studies/results/new_host_retrieval.v1/NUMERIC_QA.json)
checks all 21 pooled objectives/gradients, 131,506 selection/final scores, 968
complete-host metric rows and 486 paired contrasts. It uses independent matrix
predictions, exact-fraction no-hit calculations and expanded-host bootstrap draws.
The recorded run took 19.4 seconds with one CPU thread and about 164 MiB peak RSS;
it made no genome downloads, GenoPHI fits or paid inference calls.

Use Python 3.12.14 and the [numerical constraints](../data/studies/new_host_retrieval_environment.v1.txt).
The optional fetch obtains only the two pinned annotation tables:

```bash
pip install -e . -c data/studies/new_host_retrieval_environment.v1.txt
python scripts/run_new_host_retrieval.py \
  --fetch --annotation-data-root artifacts/annotation_inputs \
  --code-ref 6377eda61c50350a66995f6371e97080b3f026c3 \
  --out artifacts/new_host_retrieval_reproduction
python scripts/verify_new_host_retrieval.py \
  --bundle artifacts/new_host_retrieval_reproduction \
  --annotation-data-root artifacts/annotation_inputs \
  --out artifacts/new_host_retrieval_reproduction_QA.json
```

Outputs must be new. GenoPHI retains its transductive genome-vocabulary scope and
inherited environment. All scores target source MLC positivity; they do not
establish productive infection, cocktail efficacy or treatment benefit.
