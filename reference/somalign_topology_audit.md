# Compute a persistent-homology topology audit for a somalign fit

Computes the H0 (connected-component) persistence diagram of the query,
corrected-query, and reference codebooks and reports how many robustly
separated populations survive at a given distance threshold. Batch
correction that merges or erases a population shows up as a drop in the
number of H0 components between the query and corrected codebooks.

## Usage

``` r
somalign_topology_audit(
  fit,
  threshold = NULL,
  use_tda = FALSE,
  nodes = c("correction_allowed", "all")
)
```

## Arguments

- fit:

  A `somalign_fit` object.

- threshold:

  Numeric scalar in reference-scaled Euclidean units. H0 components with
  persistence (death - birth) greater than this value are counted as
  robustly separated populations. `NULL` (default) derives the threshold
  from the median of `reference$distance_quantiles`'s 95th percentile
  column – the reference's own within-population distance spread.

- use_tda:

  Logical. When `TRUE`, additionally compute H0 + H1 via
  [`TDA::ripsDiag()`](https://rdrr.io/pkg/TDA/man/ripsDiag.html) if the
  `TDA` package is installed (falls back silently to the base-R H0-only
  result, with a one-time message, when `TDA` is absent). Default
  `FALSE`.

- nodes:

  One of `"correction_allowed"` (default; recommended – nodes with no
  correction contribute unchanged, potentially spurious topology) or
  `"all"`. Note that `somalign_epsilon_sweep(..., topology = TRUE)`
  reports its topology columns using `"all"` (the correction-allowed set
  is epsilon-dependent, so it cannot be held fixed across a sweep); pass
  `nodes = "all"` here to reproduce those numbers for a given fit.

## Value

A list of class `somalign_topology`; see the source fields `threshold`,
`threshold_source` (`"auto"`/`"user"`), `n_components_query`,
`n_components_corrected`, `n_components_reference`, `topology_delta`
(corrected minus query; negative means merging),
`diagram_query`/`diagram_corrected`/`diagram_reference` (data frames of
birth/death/persistence), `bottleneck_h0` and `tda_*` (`NULL` unless
`use_tda = TRUE` and `TDA` is installed), and `topology_warning` (`TRUE`
when `topology_delta != 0`; a warning is also emitted).

## Details

This is a pure diagnostic: it reads `fit` and returns a new object
without modifying it. It is not run automatically inside
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
(too expensive to compute by default); call it directly or via
`somalign_diagnostics(fit, topology = TRUE)`.

High-dimensional marker spaces (p \> 20) are valid but nearest-neighbour
distances concentrate, which can make H0 components less visually
intuitive – a known property of topological data analysis in high
dimensions, not a bug.

## Role

This is a **diagnostic on the correction path**, answering "did the
barycentric correction merge distinct populations?" – i.e. is the
*corrected codebook* trustworthy as batch-corrected coordinates. It does
**not** bear on label transfer, which is computed from the transport
plan and never uses `node_shifts`. Treat a large negative
`topology_delta` as a reason not to use the corrected coordinates for
downstream re-clustering, not as a problem with the transferred labels.

## See also

[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md),
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
fit <- somalign_fit(qry, ref)
#> somalign_fit: 1 query node(s) have match_mass_ratio > 1 (max 1.09); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_topology_audit(fit)
#> <somalign_topology>
#>   threshold: 0.8892 (auto)
#>   H0 components  query: 3  corrected: 3  reference: 4
#>   topology_delta: +0   warning: FALSE
```
