# Tune transport-plan knobs against label-transfer accuracy

The supervised counterpart to
[`somalign_select_epsilon()`](https://mdmanurung.github.io/somalign/reference/somalign_select_epsilon.md):
instead of an unsupervised plan-geometry criterion, it selects
transport-plan parameters by cross-validated label-transfer performance.
For each parameter combination it runs stratified k-fold CV (reusing
pre-trained SOMs per fold for efficiency – plan knobs change only the OT
solve, not the SOMs) and scores pooled held-out predictions with
[`somalign_label_metrics()`](https://mdmanurung.github.io/somalign/reference/somalign_label_metrics.md).

## Usage

``` r
somalign_tune(
  data,
  labels,
  grid,
  param_grid,
  k = 5L,
  metric = c("mcc", "macro_f1", "accuracy", "ece"),
  stratify = TRUE,
  rlen = 20L,
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  solver = "internal",
  max_iter = 1000,
  tol = 1e-07,
  n_bins = 10L,
  seed = 1L
)
```

## Arguments

- data:

  Numeric cell-by-feature matrix.

- labels:

  Character vector of per-cell labels.

- grid:

  A
  [`kohonen::somgrid`](https://rdrr.io/pkg/kohonen/man/unit.distances.html).

- param_grid:

  A data frame (one row per combination of the scalar knobs `epsilon`,
  `rho_query`, `rho_ref`, `diagonal_boost`) or a list of named lists
  (which may additionally carry `feature_weights`). Each must specify
  `epsilon`.

- k:

  Folds. Default `5L`.

- metric:

  Objective to optimise: `"mcc"` (default), `"macro_f1"`, `"accuracy"`
  (all maximised) or `"ece"` (minimised).

- stratify, rlen, seed:

  As in
  [`somalign_cross_validate()`](https://mdmanurung.github.io/somalign/reference/somalign_cross_validate.md).

- min_match_fraction, confidence_threshold:

  Label-acceptance gates, matching
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  defaults.

- solver, max_iter, tol:

  OT solver settings.

- n_bins:

  Calibration bins for the `ece` column. Default `10L`.

## Value

A list of class `somalign_tune` with `best` (the winning combo as a
one-row data frame), `best_params` (named list), `grid` (all
combinations with their CV metrics), and `metric`.

## Details

Tunable knobs are those that shape the transport plan without needing
anchors or query labels: `epsilon`, `rho_query`, `rho_ref`,
`diagonal_boost`, and `feature_weights` (a numeric per-feature vector).
`label_guided` and `rho_anchor` are out of scope here (they require a
labelled query SOM or anchor pairs, respectively).

Note that `"mcc"` and `"accuracy"` are scored on *accepted* predictions
only, so they can be inflated by settings that abstain on hard cells
(higher `epsilon` raises accuracy while dropping `coverage`). Always
read the `coverage` and `macro_f1` columns alongside the objective:
`macro_f1` falls when rare classes are abstained away, making it a more
coverage-robust target for imbalanced data.

## Examples

``` r
# \donttest{
if (requireNamespace("kohonen", quietly = TRUE)) {
  set.seed(1)
  x <- rbind(matrix(rnorm(200, -2), ncol = 2), matrix(rnorm(200, 2), ncol = 2))
  colnames(x) <- c("f1", "f2")
  lab <- rep(c("low", "high"), each = 100)
  tuned <- somalign_tune(x, lab, grid = kohonen::somgrid(2, 2, "hexagonal"),
    param_grid = data.frame(epsilon = c(0.05, 0.1, 0.2)), k = 3)
  tuned$best_params$epsilon
}
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> [1] 0.05
# }
```
