# Return worst-projecting query SOM nodes

Returns the `n` query nodes with the lowest match fraction — the nodes
whose mass the OT solver could not route to the reference — sorted
ascending. Includes the dominant reference label each node maps to.

## Usage

``` r
somalign_worst_nodes(fit, n = 10L)
```

## Arguments

- fit:

  A `somalign_fit` object.

- n:

  Number of nodes to return. Default `10`.

## Value

A data frame with columns `query_node`, `query_mass`,
`transported_mass`, `match_fraction`, `correction_allowed`,
`correction_norm`, and `top_ref_label`, one row per node.

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
somalign_worst_nodes(fit, n = 4)
#>    query_node query_mass transported_mass match_fraction correction_allowed
#> V2          2        0.3        0.2662839      0.8876129               TRUE
#> V3          3        0.5        0.4460688      0.8921376               TRUE
#> V1          1        0.1        0.1179511      1.0000000               TRUE
#> V4          4        0.1        0.1069541      1.0000000               TRUE
#>    correction_norm top_ref_label
#> V2       0.1056234          <NA>
#> V3       0.4390426          <NA>
#> V1       0.3455129          <NA>
#> V4       0.1478819          <NA>
```
