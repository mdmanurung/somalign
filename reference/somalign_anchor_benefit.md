# Quantify the label-transfer benefit of anchor (repeat) samples

Measures how much anchor samples – QC/repeat specimens run in both the
reference and query batches – improve label transfer, by sweeping the
anchor cost-bonus strength `rho_anchor` and scoring the transferred
labels against a *known* query label. `rho_anchor = 0` is the exact
no-anchor baseline (plain
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md));
positive values bias the transport plan toward anchor-supported node
pairs, exactly as `somalign_fit_anchored(correction = "cost_bonus")`.
Anchors influence labels only through this cost bonus, so this single
sweep captures their entire label-transfer effect (the `"subspace"`
correction leaves labels identical to `rho_anchor = 0`).

## Usage

``` r
somalign_anchor_benefit(
  query,
  reference,
  query_labels,
  anchor_old,
  anchor_new,
  rho_grid = c(0, 0.5, 1, 2, 5),
  metric = c("mcc", "macro_f1", "accuracy"),
  eval_mask = NULL,
  epsilon = 0.1,
  rho_query = 1,
  rho_ref = 1,
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  solver = "internal",
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L,
  n_bins = 10L
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A labelled `somalign_reference` object.

- query_labels:

  Character ground-truth labels, one per query cell.

- anchor_old, anchor_new:

  Paired anchor cell matrices (repeat samples in the reference and query
  batches respectively), same number of rows.

- rho_grid:

  Non-negative anchor strengths to sweep. Default `c(0, 0.5, 1, 2, 5)`;
  include `0` for the no-anchor baseline.

- metric:

  Objective for `best`/`lift`: `"mcc"` (default), `"macro_f1"`, or
  `"accuracy"` (all maximised).

- eval_mask:

  Optional logical vector, one per query cell, selecting the cells to
  score (e.g. exclude the anchor samples for a clean held-out measure).
  Default `NULL` scores all cells.

- epsilon, rho_query, rho_ref, solver, max_iter, tol:

  OT settings.

- min_match_fraction, confidence_threshold:

  Label-acceptance gates.

- chunk_size:

  Anchor-projection chunk size.

- n_bins:

  Calibration bins for the `ece` column.

## Value

A list of class `somalign_anchor_benefit` with `grid` (one row per
`rho_anchor`: accuracy, macro_f1, mcc, coverage, ece), `baseline` (the
`rho_anchor = 0` row), `best`, `lift` (best minus baseline on `metric`),
and `metric`.

## Details

This is a *validation* tool: it needs ground-truth query labels (e.g. an
independently gated held-out batch), which production label transfer
does not have. It reuses the fixed reference/query SOMs and a single
anchor count matrix across the sweep, recomputing only the OT solve, so
it is cheap.

## See also

[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md),
[`somalign_cross_validate()`](https://mdmanurung.github.io/somalign/reference/somalign_cross_validate.md)
