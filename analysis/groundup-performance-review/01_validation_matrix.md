# Validation Matrix For The Next Improvement Round

Date: 2026-07-19

This matrix turns the multi-team discussion into falsifiable experiments. The
first seven rows are the package-performance approaches to validate alone and in
combination; the eighth row is engineering performance.

## Shared Reporting Contract

Every experiment should write:

- `summary.csv`: one row per method/seed/split with all primary metrics.
- `per_class_metrics.csv`: precision, recall, F1, support, and abstention-aware
  counts.
- `predictions.rds` or a compact equivalent: per-cell truth, prediction,
  confidence, margin, transferred posterior, source method, sample ID, and fold.
- `frequency_correlations.csv`: hard and soft repeat-sample CLR weighted-r where
  sample IDs are available.
- `calibration.csv`: ECE, MCE, Brier, reliability bins, and class-stratified
  calibration.
- `provenance.json`: package version, git SHA, command, seed, input files,
  artifact hashes, R session, host, and timestamp.

Primary label-transfer metrics must be `accuracy_all`, full macro-F1, MCC,
coverage, per-class precision/recall/F1, and confusion with abstain counted.
Accepted-cell metrics can be secondary only.

## Experiment 1: Multi-Objective Epsilon/Rho Selection

Hypothesis: a predeclared utility that penalizes abstention and poor calibration
selects settings that generalize better than accepted-cell MCC.

Design:

- Sweep `epsilon = c(0.02, 0.05, 0.1, 0.2)`.
- Sweep practical `rho_query`/`rho_ref` values around the current defaults.
- Evaluate held-out folds, Batch1-to-Batch2, and repeat-sample frequencies.
- Compare selection objectives: accepted MCC, full macro-F1, calibration-aware
  utility, and OtherT-constrained utility.

Pass:

- Full macro-F1 and `accuracy_all` do not regress.
- Coverage remains within a predeclared floor.
- ECE/Brier improve or remain stable.
- OtherT recall improves without precision dropping below the fixed threshold.

## Experiment 2: Calibrated Acceptance

Hypothesis: label posterior, margin, entropy, match fraction, and transported
mass can produce calibrated accept/rescue decisions.

Design:

- Fit calibration models on reference-only CV folds.
- Test isotonic, logistic, and temperature-style calibrators.
- Report reliability by class and by confidence bin.
- Evaluate fixed precision curves and coverage curves.

Pass:

- ECE and Brier improve on held-out folds.
- Fixed precision thresholds transfer to Batch1-to-Batch2 without retuning.
- Abstention reduction does not degrade full macro-F1.

## Experiment 3: Strict T/NK Router Plus Local Specialist Projection

Hypothesis: OtherT is a T/NK boundary problem and should be handled by a
high-specificity local model.

Design:

- Route only deployable T/NK candidates using marker evidence, global labels,
  margins, entropy, soft probabilities, and outside-reference diagnostics.
- Exclude B/myeloid/plasmacytoid evidence explicitly.
- Run local T/NK projection with T/NK-focused features.
- Test the router alone, local projection alone on routed cells, and the
  combined replacement rule.

Pass:

- Non-T/NK false-route rate is sharply reduced versus the broad router.
- B-cell and Mono/DC recall do not collapse.
- OtherT recall increases at fixed precision.
- Full macro-F1 improves or remains stable.

## Experiment 4: High-Margin OtherT Rescue

Hypothesis: a conservative rescue layer can recover true OtherT cells missed by
the high-precision baseline.

Design:

- Keep `epsilon = 0.1` baseline predictions as default.
- Generate candidates from `epsilon = 0.02` and T-marker-weighted
  `epsilon = 0.05`.
- Rescue only cells with strong OtherT posterior, T-cell evidence, ambiguous
  baseline margin/abstention, and no clean CD4/CD8/NK signature.
- Sweep posterior and margin thresholds.

Pass:

- OtherT recall increases while precision remains >= baseline tolerance or a
  predeclared fixed threshold.
- Full macro-F1 does not decline.
- Improvements reproduce at 120k scale and across seeds.

## Experiment 5: Shared Soft Posterior

Hypothesis: using per-cell soft posteriors for both labels and frequencies
reduces node-boundary artifacts.

Design:

- Compute `cell -> query-node softness -> OT correspondence -> reference label
  probabilities`.
- Compare pure OT posterior, direct reference-kNN posterior, and blended
  posterior.
- Use the same posterior for hard labels, abstention, and sample frequencies.

Pass:

- Repeat-sample metacluster CLR weighted-r equals or exceeds current 0.9566.
- Hard-label macro-F1 and OtherT F1 do not regress.
- Calibration improves versus node-argmax labels.

## Experiment 6: Reference-Only Feature/Metric Learning

Hypothesis: data-driven feature weights can improve difficult class boundaries
without using query truth.

Design:

- Learn diagonal weights from reference CV using class separation,
  within-class compactness, rare-class preservation, and shrinkage.
- Compare unweighted, manual T-marker weights, anchor-derived weights, and
  learned weights.
- Test stability across folds, seeds, and bootstrap samples.

Pass:

- Learned weights beat both unweighted and manual T-marker baselines on
  held-out folds.
- OtherT F1 improves without major loss in non-T cell classes.
- Weight profiles are biologically interpretable and stable.

## Experiment 7: Refined OtherT Taxonomy

Hypothesis: OtherT is too heterogeneous as one residual label; transferring
substates and collapsing back improves boundary behavior.

Design:

- Subcluster reference OtherT with T/NK markers only.
- Annotate candidate sublabels by marker profiles.
- Transfer sublabels locally, then collapse to OtherT.
- Include a negative-control random sublabel split.

Pass:

- Collapsed OtherT precision and recall improve.
- Sublabels have stable marker profiles.
- Random splits do not improve performance.

## Experiment 8: Runtime/API Performance

Hypothesis: projection memory and object retention can be reduced without
changing scientific outputs.

Design:

- Add and test constructor-level chunking.
- Reuse `som$unit.classif` where valid.
- Stream corrected projection instead of materializing corrected matrices.
- Add compact storage controls.
- Benchmark `somalign_query()`, `somalign_reference()`, `somalign_fit()`, and
  `somalign_results()` across cells, markers, nodes, and chunk sizes.

Pass:

- Outputs are equivalent under deterministic tests.
- Peak memory drops on large inputs.
- Runtime does not regress materially for small inputs.

## Combination Tests

After single-method validation, run these combinations:

1. Multi-objective tuning + calibrated acceptance.
2. Strict T/NK router + local specialist projection + high-margin rescue.
3. Shared soft posterior + calibrated acceptance.
4. Feature-weight learning + strict T/NK router.
5. Refined OtherT taxonomy + local specialist projection.
6. Best label-transfer method + tuned soft-frequency estimator.
7. Best scientific method + runtime/API hardening equivalence tests.

## Stop Conditions

Reject any method that:

- improves accepted-cell accuracy while reducing full macro-F1 or coverage,
- improves OtherT recall only by collapsing non-T classes,
- uses query truth or truth-derived labels in a deployable path,
- improves repeat correlation only by suppressing rare clusters,
- lacks reproducible artifacts and provenance,
- changes runtime behavior without equivalence tests.
