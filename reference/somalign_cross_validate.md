# Cross-validate label transfer on held-out cells

Stratified k-fold cross-validation that measures how accurately somalign
transfers labels to cells it has not seen. Each fold trains a reference
SOM on the training split, projects the held-out split, and transfers
labels; every held-out cell has a real ground-truth label, so this
estimates true generalisation without any external labelled query data.
Pooled results are scored with
[`somalign_label_metrics()`](https://mdmanurung.github.io/somalign/reference/somalign_label_metrics.md)
and
[`somalign_calibration()`](https://mdmanurung.github.io/somalign/reference/somalign_calibration.md).

## Usage

``` r
somalign_cross_validate(
  data,
  labels,
  grid,
  k = 5L,
  stratify = TRUE,
  epsilon = 0.1,
  solver = "internal",
  rlen = 20L,
  n_bins = 10L,
  seed = 1L,
  ...
)
```

## Arguments

- data:

  Numeric cell-by-feature matrix.

- labels:

  Character vector of per-cell labels, `nrow(data)` long.

- grid:

  A
  [`kohonen::somgrid`](https://rdrr.io/pkg/kohonen/man/unit.distances.html)
  for the reference and query SOMs.

- k:

  Positive integer. Number of folds. Default `5L`.

- stratify:

  Logical. Stratify folds by class. Default `TRUE`.

- epsilon, solver:

  Passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- rlen:

  SOM training iterations for both SOMs. Default `20L`.

- n_bins:

  Calibration bins. Default `10L`.

- seed:

  Integer or `NULL`. RNG seed, restored on exit. Default `1L`.

- ...:

  Further arguments forwarded to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## Value

A list of class `somalign_cross_validation` with `metrics`
([`somalign_label_metrics()`](https://mdmanurung.github.io/somalign/reference/somalign_label_metrics.md)),
`calibration`
([`somalign_calibration()`](https://mdmanurung.github.io/somalign/reference/somalign_calibration.md)),
`per_fold` (data frame), `predictions` (pooled per-cell data frame),
`k`.

## Examples

``` r
# \donttest{
if (requireNamespace("kohonen", quietly = TRUE)) {
  set.seed(1)
  x <- rbind(matrix(rnorm(200, -2), ncol = 2), matrix(rnorm(200, 2), ncol = 2))
  colnames(x) <- c("f1", "f2")
  lab <- rep(c("low", "high"), each = 100)
  cv <- somalign_cross_validate(x, lab,
    grid = kohonen::somgrid(2, 2, "hexagonal"), k = 3)
  cv$metrics$accuracy
}
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.06); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.11); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.27); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> [1] 0.985
# }
```
