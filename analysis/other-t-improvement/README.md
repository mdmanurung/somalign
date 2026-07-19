# OtherT Label-Transfer Improvement Validation

This folder validates seven publication-readiness improvement ideas for
somalign's BMV pilot label-transfer use case, with emphasis on the `OtherT`
compartment and repeat-sample SOM-code frequency correlations.

The analysis is intentionally outside the package build. The root
`.Rbuildignore` excludes `analysis/`, so these scripts document and validate
claims without becoming installed package content.

## Run

Fast pilot:

```sh
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  analysis/other-t-improvement/other_t_improvement_experiment.R 30000 10 20
```

Notebook true-repeat audit:

```sh
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  analysis/other-t-improvement/notebook_repeat_sample_audit.R
```

Fuller 120k-style run:

```sh
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  analysis/other-t-improvement/other_t_improvement_experiment.R 120000 15 40
```

Arguments are positional:

1. `n_subsample`: total stratified cells across reference/query batches.
2. `grid_side`: uses a square `grid_side x grid_side` hexagonal SOM.
3. `rlen`: SOM training iterations.
4. `seed`: optional, defaults to `1`.

Environment overrides:

- `SOMALIGN_PILOT_QS`: labelled BMV pilot `.qs` file.
- `SOMALIGN_REF_BATCH`: reference batch, default `Batch1`.
- `SOMALIGN_QRY_BATCH`: labelled query batch, default `Batch2`.
- `SOMALIGN_LABEL_COL`: label column, default `gate_lineage2`.
- `SOMALIGN_BATCH_COL`: batch column, default `batch_id`.
- `SOMALIGN_NOTEBOOK_PROJECT_DIR`: external example-notebook project, default
  `/exports/para-lipg-hpc/mdmanurung/bmv_pilot_cytof_integration`.
- `SOMALIGN_NOTEBOOK_SOMALIGN_DIR`: external notebook SOMalign artifact folder.
- `SOMALIGN_REPEAT_CSV`: true repeat metadata table, default
  `repeated_samples.csv` in the notebook project.

## What Is Validated

1. Soft abundance frequencies versus hard nearest-node frequencies for repeated
   sample or pseudo-repeat SOM-code and lineage profiles.
2. A targeted `OtherT` validation grid scored with full per-class precision,
   recall and F1, counting abstentions as false negatives.
3. The `epsilon = 0.2` operating point against the shipped `epsilon = 0.1`.
4. Marker-aware transport costs, using explicit T/NK boundary marker weights
   and anchor-derived feature weights when anchors are available.
5. A T/NK-local projection stage and a production-like local rerouting combo.
6. Margin and transport-entropy triage, reported separately from recall because
   selective acceptance usually trades recall for precision.
7. Anchor strength tuning over `rho_anchor`, including default `rho_anchor = 1`
   and best observed settings.

## Outputs

- `run_metadata.csv`: dataset, batch, grid and repeat-group metadata.
- `label_summary.csv`: aggregate label-transfer metrics per approach/config.
- `per_class_metrics.csv`: full per-class precision, recall and F1.
- `abundance_correlation_pairs.csv`: per-repeat hard/soft CLR correlations.
- `abundance_correlation_summary.csv`: mean/median/weighted hard/soft
  correlation summary.
- `notebook_repeat_sample_metadata.csv`: the 20 true repeat sample IDs from the
  SOMalign example notebook, with pilot and new/BMV FCS filenames.
- `notebook_repeat_correlation_pairs.csv`: per-sample hard/soft CLR weighted-r
  from the notebook true-repeat artifact.
- `notebook_repeat_correlation_summary.csv`: hard versus soft median CLR
  weighted-r at metacluster and lineage resolution for true repeats.
- `notebook_repeat_soft_grid_summary.csv`: k-grid summary from the dedicated
  soft-frequency repeat-sample artifact.
- `notebook_repeat_source_artifacts.csv`: source paths used by the audit.
- `targeted_grid_ranked.csv`: deployable global parameter rows ranked for
  `OtherT`.
- `targeted_grid_ranked_all.csv`: deployable rows plus oracle upper-bound rows.
- `router_diagnostics.csv`: T/NK reroute precision/recall diagnostics.
- `anchor_grid.csv`: `rho_anchor` sweep for proxy same-lineage anchors and a
  random mispaired negative control.
- `other_t_improvement_results.rds`: all objects needed for reinspection.
- `session_info.txt`: R session provenance.

## Interpretation Rules

Use `evidence_class` before making claims:

- `deployable`: can be implemented without query truth.
- `oracle_upper_bound`: uses held-out query labels to build query SOM label
  probabilities; useful only as an upper bound.
- `oracle_subset`: selects true T/NK query cells before local projection; useful
  for compartment-capacity testing, not a production workflow.
- `proxy_oracle_anchor`: pairs anchor cells by true lineage because true repeat
  metadata were unavailable; useful as a proxy stress test.
- `pseudo_repeat_control`: balanced pseudo-repeat groups, not true repeated
  samples. Do not use these rows to support publication claims about repeat
  sample reproducibility.
- `true_repeat_metadata`: the notebook's pilot-A versus BMV/query-B repeat
  design joined by `repeated_samples.csv`; use these rows for repeat-frequency
  publication claims.

For publication claims, headline `accuracy_all`, `full_macro_f1`,
`otherT_precision`, `otherT_recall`, `otherT_f1`, and `coverage`. Accepted-cell
metrics are secondary because abstentions can inflate them.
