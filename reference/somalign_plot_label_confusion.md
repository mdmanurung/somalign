# Plot label transfer confusion heatmap

Row-normalised heatmap of old-to-transferred label pairs for accepted
cells. High values on the diagonal indicate coherent transfer; strong
off-diagonal entries warrant further inspection.

## Usage

``` r
somalign_plot_label_confusion(fit, min_confidence = NULL)
```

## Arguments

- fit:

  A `somalign_fit` object.

- min_confidence:

  Minimum `transferred_label_confidence` to include. `NULL` (default)
  imposes no additional filter beyond acceptance.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
fit <- somalign_fit(qry, ref)
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.18); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_plot_label_confusion(fit)
#> Error: No accepted transferred labels found; cannot build confusion plot.
```
