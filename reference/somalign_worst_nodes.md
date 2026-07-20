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
fit <- somalign_fit(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_worst_nodes(fit, n = 4)
#>    query_node query_mass transported_mass match_fraction correction_allowed
#> V2          2        0.3        0.2869031      0.9563437               TRUE
#> V1          1        0.1        0.1234245      1.0000000               TRUE
#> V3          3        0.5        0.5002755      1.0000000               TRUE
#> V4          4        0.1        0.1111685      1.0000000               TRUE
#>    correction_norm transport_entropy top_ref_label
#> V2       0.1185567      6.652620e-01          <NA>
#> V1       0.3455139      3.089021e-02          <NA>
#> V3       0.4395087      9.990417e-01          <NA>
#> V4       0.1478769      8.307340e-10          <NA>
```
