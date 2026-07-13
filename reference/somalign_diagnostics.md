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
fit <- somalign_fit(qry, ref)
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.91); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
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
#> [1] 20
#> 
#> $solver$converged
#> [1] TRUE
#> 
#> $solver$final_delta
#> [1] 7.490753e-08
#> 
#> $solver$epsilon
#> [1] 0.5
#> 
#> $solver$rho_query
#> [1] 1
#> 
#> $solver$rho_ref
#> [1] 1
#> 
#> $solver$cost_scale
#> [1] 1.660105
#> 
#> $solver$rel_marginal_row_error
#> [1] 0.09492421
#> 
#> $solver$rel_marginal_col_error
#> [1] 0.08788304
#> 
#> 
#> $ot
#> $ot$transport_mass
#> [1] 1.288971
#> 
#> $ot$row_mass
#>        V1        V2        V3        V4 
#> 0.1914436 0.3949242 0.5463233 0.1562799 
#> 
#> $ot$col_mass
#>        V1        V2        V3        V4 
#> 0.3756875 0.2878830 0.4683774 0.1570230 
#> 
#> $ot$query_mass
#> [1] 0.1 0.3 0.5 0.1
#> 
#> $ot$reference_mass
#> [1] 0.3 0.2 0.4 0.1
#> 
#> $ot$match_fraction
#> [1] 1 1 1 1
#> 
#> $ot$match_mass_ratio
#> [1] 1.914436 1.316414 1.092647 1.562799
#> 
#> $ot$max_row_mass_error
#> [1] 0.09492421
#> 
#> $ot$max_col_mass_error
#> [1] 0.08788304
#> 
#> 
#> $nodes
#>    query_node query_mass transported_mass match_fraction correction_allowed
#> V1          1        0.1        0.1914436              1               TRUE
#> V2          2        0.3        0.3949242              1               TRUE
#> V3          3        0.5        0.5463233              1               TRUE
#> V4          4        0.1        0.1562799              1               TRUE
#>    correction_norm
#> V1       0.4361949
#> V2       0.4044215
#> V3       0.5557118
#> V4       0.4532371
#> 
#> $projection
#> $projection$outside_direct_fraction
#> [1] 0.3
#> 
#> $projection$outside_corrected_fraction
#> [1] 0.3
#> 
#> 
```
