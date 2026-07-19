# OtherT Improvement Validation Findings

Run: `other_t_improvement_experiment.R 30000 10 20`

Data: labelled BMV pilot, `Batch1` as reference and `Batch2` as labelled query.
The run used 30,000 cells balanced across batch and `gate_lineage2` labels,
a 10 x 10 hexagonal SOM, `rlen = 20`, and seed 1.

## Repeat-Sample Frequency Correlation

The true repeat metadata are in the SOMalign example notebook project at
`/exports/para-lipg-hpc/mdmanurung/bmv_pilot_cytof_integration/repeated_samples.csv`.
That table contains 20 matched sample IDs, each with one pilot acquisition and
one new/BMV acquisition. The relevant design is therefore pilot measurement A
versus BMV/query projection B, not the labelled `Batch1` versus `Batch2` split
used for the OtherT label-transfer pilot.

The notebook artifacts validate the repeat-frequency claim:

| resolution | hard nearest-node median CLR weighted r | soft kNN median CLR weighted r | median difference delta | samples soft > hard | evidence |
|---|---:|---:|---:|---:|---|
| metacluster | 0.9322 | 0.9566 | +0.0243 | 20/20 | true repeat metadata |
| lineage | 0.9873 | 0.9938 | +0.0065 | 14/20 | true repeat metadata |

The dedicated k-grid artifact is consistent: `soft k=8` gives median
metacluster CLR weighted r 0.9564 versus 0.9322 for hard nearest-node, and
soft beats hard in 20/20 repeat samples. The split-half sampling ceiling is
0.9985, so the soft kNN gain closes a meaningful fraction of the remaining
metacluster headroom. The median paired per-sample metacluster delta is +0.0192.
This is now suitable evidence for the package claim, provided the source
artifact paths are retained or the audit is rerun before release.

The labelled Batch1/Batch2 pilot still does not contain this true repeat
metadata. Its pseudo-repeat control remains a negative/control artifact:

| scope | hard direct weighted r | soft kNN weighted r | evidence |
|---|---:|---:|---|
| SOM unit | 0.8719 | 0.5930 | pseudo-repeat control |
| lineage | NA | NA | no finite pairs |

Do not use the pseudo-repeat rows for publication claims about repeat-sample
reproducibility.

## Deployable OtherT Results

Baseline (`epsilon = 0.1`) gives high OtherT precision but low recall:

| config | coverage | accuracy_all | full macro-F1 | OtherT precision | OtherT recall | OtherT F1 |
|---|---:|---:|---:|---:|---:|---:|
| global_epsilon_0.1 | 0.9462 | 0.9215 | 0.9441 | 0.9653 | 0.7236 | 0.8272 |

Best deployable settings improved OtherT recall and F1, but not precision:

| config | coverage | accuracy_all | full macro-F1 | OtherT precision | OtherT recall | OtherT F1 |
|---|---:|---:|---:|---:|---:|---:|
| global_epsilon_0.02 | 0.9923 | 0.9563 | 0.9596 | 0.9140 | 0.8632 | 0.8879 |
| t_marker_weights_epsilon_0.05 | 0.9874 | 0.9545 | 0.9600 | 0.9306 | 0.8528 | 0.8900 |
| baseline_plus_global_epsilon_0.02_otherT_margin_ge_0.00 | 0.9792 | 0.9451 | 0.9544 | 0.9135 | 0.8656 | 0.8889 |
| baseline_plus_t_marker_weights_epsilon_0.05_otherT_margin_ge_0.00 | 0.9740 | 0.9430 | 0.9546 | 0.9306 | 0.8528 | 0.8900 |
| global_epsilon_0.05_rho_2_2 | 0.9779 | 0.9487 | 0.9583 | 0.9542 | 0.8176 | 0.8807 |

No deployable configuration or rescue combination in this pilot improved both
OtherT precision and recall over the baseline simultaneously. The best
practical choice depends on the operating goal:

- Prefer `t_marker_weights_epsilon_0.05` for the best OtherT F1 and macro-F1
  while retaining moderate precision.
- Prefer `global_epsilon_0.02` when recall and coverage matter more than
  precision.
- Keep `epsilon = 0.1` if high OtherT precision is the priority.

## Approach Status

1. Soft frequencies: validated for true repeats in the SOMalign example
   notebook artifacts. At metacluster resolution, hard nearest-node median CLR
   weighted r was 0.9322 and soft kNN was 0.9566, with soft winning 20/20
   samples.
2. Targeted OtherT grid: useful. Low epsilon and T-marker weights gave the best
   deployable F1/recall improvements.
3. `epsilon = 0.2`: not useful for OtherT recall. It raised precision to
   0.9994 but dropped recall to 0.6568.
4. Marker-aware transport cost: useful. Explicit T/NK marker weights at
   `epsilon = 0.05` were the best deployable F1 result.
5. T/NK-local projection: promising only with a better router. Truth-selected
   local projection improved compartment capacity, but the deployable router was
   too broad: it routed 97.7% of query cells and 93.4% of non-T/NK cells.
6. Margin/entropy triage: useful as QC, not as a recall-improving strategy.
   Margin thresholds did not change the baseline in this pilot; entropy gating
   reduced coverage and full recall.
7. Anchor tuning: default `rho_anchor = 1` was inert. Proxy same-lineage anchors
   and random mispaired controls showed no label-transfer lift in this real
   Batch1-to-Batch2 pilot.

The deployable rescue combinations, which keep the baseline and selectively
rescue OtherT calls from `global_epsilon_0.02` or
`t_marker_weights_epsilon_0.05`, improved recall but did not preserve baseline
precision.

## Recommended Next Experiment

The most plausible route to improving both precision and recall is a stricter
T/NK router followed by local projection, not a global transport-only knob.
The current router uses global top or second label and is far too broad. Next,
test a router based on high-specificity direct labels plus outside-reference,
marker-gate, and margin conditions, then reroute only the ambiguous T/NK
boundary cells.
