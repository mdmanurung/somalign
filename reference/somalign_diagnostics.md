# Extract somalign diagnostics

Extract somalign diagnostics

## Usage

``` r
somalign_diagnostics(fit)
```

## Arguments

- fit:

  A `somalign_fit` object.

## Value

A named list of solver, OT, node, and projection diagnostics.

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
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_diagnostics(fit)
#> $solver
#> $solver$requested
#> [1] "internal"
#> 
#> $solver$used
#> [1] "internal"
#> 
#> $solver$notes
#> character(0)
#> 
#> $solver$iterations
#> [1] 72
#> 
#> $solver$converged
#> [1] TRUE
#> 
#> $solver$final_delta
#> [1] 9.290696e-08
#> 
#> $solver$epsilon
#> [1] 0.1
#> 
#> $solver$rho_query
#> [1] 1
#> 
#> $solver$rho_ref
#> [1] 1
#> 
#> $solver$cost_scale
#> [1] 2.762185
#> 
#> $solver$rel_marginal_row_error
#> [1] 0.02342453
#> 
#> $solver$rel_marginal_col_error
#> [1] 0.03828297
#> 
#> 
#> $ot
#> $ot$transport_mass
#> [1] 1.021772
#> 
#> $ot$row_mass
#>        V1        V2        V3        V4 
#> 0.1234245 0.2869031 0.5002755 0.1111685 
#> 
#> $ot$col_mass
#>        V1        V2        V3        V4 
#> 0.2894307 0.2382830 0.3828895 0.1111685 
#> 
#> $ot$query_mass
#> [1] 0.1 0.3 0.5 0.1
#> 
#> $ot$reference_mass
#> [1] 0.3 0.2 0.4 0.1
#> 
#> $ot$match_fraction
#> [1] 1.0000000 0.9563437 1.0000000 1.0000000
#> 
#> $ot$match_mass_ratio
#> [1] 1.2342453 0.9563437 1.0005511 1.1116850
#> 
#> $ot$max_row_mass_error
#> [1] 0.02342453
#> 
#> $ot$max_col_mass_error
#> [1] 0.03828297
#> 
#> 
#> $nodes
#>    query_node query_mass transported_mass match_fraction correction_allowed
#> V1          1        0.1        0.1234245      1.0000000               TRUE
#> V2          2        0.3        0.2869031      0.9563437               TRUE
#> V3          3        0.5        0.5002755      1.0000000               TRUE
#> V4          4        0.1        0.1111685      1.0000000               TRUE
#>    correction_norm
#> V1       0.3455139
#> V2       0.1185567
#> V3       0.4395087
#> V4       0.1478769
#> 
#> $projection
#> $projection$outside_direct_fraction
#> [1] 0.3
#> 
#> $projection$outside_corrected_fraction
#> [1] 0.1
#> 
#> 
```
