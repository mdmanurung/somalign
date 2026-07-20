# Select epsilon from an epsilon-sweep curve

Sweeps a grid of epsilon values via
[`somalign_epsilon_sweep()`](https://mdmanurung.github.io/somalign/reference/somalign_epsilon_sweep.md)
and selects a recommended value using one of three criteria: the
susceptibility (critical-epsilon) peak, the mutual-information-vs-cost
elbow, or the smallest epsilon retaining a target fraction of maximum
mutual information.

## Usage

``` r
somalign_select_epsilon(
  query,
  reference,
  epsilon = c(0.02, 0.05, 0.1, 0.2, 0.5, 1),
  rho_query = 1,
  rho_ref = 1,
  method = c("critical", "elbow", "entropy_fraction"),
  entropy_fraction = 0.9,
  solver = c("log_domain", "internal", "auto", "annealing"),
  max_iter = 1000,
  tol = 1e-07,
  diagonal_boost = 0,
  label_guided = FALSE,
  parallel = FALSE,
  anneal_start = 10,
  anneal_factor = NULL,
  anneal_stages = 10L
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Numeric vector of candidate epsilon values, or `NULL` to use the
  default grid. Default `c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0)`.

- rho_query, rho_ref:

  Marginal relaxations passed to the OT solver.

- method:

  Character. `"critical"` (default) uses the susceptibility peak's
  recommended epsilon (`0.3 * epsilon_c`). `"elbow"` uses the
  max-second-difference rule on the mutual-information-vs-cost curve.
  `"entropy_fraction"` chooses the largest (most-regularized) epsilon
  that still retains `entropy_fraction * max(mutual_information)`.

- entropy_fraction:

  Numeric in (0, 1\]. Target fraction of maximum mutual information when
  `method = "entropy_fraction"`. Default `0.90`.

- solver, max_iter, tol:

  Sinkhorn solver parameters.

- diagonal_boost, label_guided, parallel:

  See
  [`somalign_epsilon_sweep()`](https://mdmanurung.github.io/somalign/reference/somalign_epsilon_sweep.md).

- anneal_start, anneal_factor, anneal_stages:

  Annealing-schedule tuning parameters, used only when
  `solver = "annealing"`. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## Value

A list of class `"somalign_epsilon_selection"` with `selected_epsilon`,
`curve` (the sweep's `table`), and `method`.

## See also

[`somalign_epsilon_sweep()`](https://mdmanurung.github.io/somalign/reference/somalign_epsilon_sweep.md)

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
somalign_select_epsilon(qry, ref, epsilon = c(0.05, 0.1, 0.2))
#> $selected_epsilon
#> [1] 0.03
#> 
#> $curve
#>   epsilon log_epsilon       Phi        log_Z mutual_information
#> 1    0.05   -2.995732 0.3382930 -0.003142188           1.248706
#> 2    0.10   -2.302585 0.3503857 -0.138063739           1.204442
#> 3    0.20   -1.609438 0.3827881 -0.410877619           1.102695
#>   conditional_entropy_mean expected_cost transport_mass cost_scale iterations
#> 1                0.3746289     0.1096809      0.9818533   2.762185        113
#> 2                0.4237985     0.1147550      1.0217716   2.762185         66
#> 3                0.5551782     0.1281361      1.1035035   2.762185         39
#>   converged susceptibility
#> 1      TRUE             NA
#> 2      TRUE      0.0320964
#> 3      TRUE             NA
#> 
#> $method
#> [1] "critical"
#> 
#> $sweep
#> somalign epsilon sweep [3 points]
#>   epsilon range      : 0.05 - 0.2
#>   critical epsilon   : 0.1
#>   recommended epsilon (0.3x critical): 0.03
#>   cost_scale         : 2.762
#> 
#> attr(,"class")
#> [1] "somalign_epsilon_selection"
```
