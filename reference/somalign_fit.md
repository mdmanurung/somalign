# Align a query SOM to a reference SOM

Align a query SOM to a reference SOM

## Usage

``` r
somalign_fit(
  query,
  reference,
  epsilon = 0.5,
  rho_query = 1,
  rho_ref = 1,
  solver = c("internal", "log_domain", "auto"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Entropic regularisation strength. The cost matrix is normalised by its
  median positive entry before computing the Sinkhorn kernel, so
  `epsilon` is approximately scale- and dimension-invariant. Values
  around `0.5` give meaningful regularisation for typical z-scored SOM
  codebooks; very small values (\< 0.1) make the transport increasingly
  discrete. The normalisation scale is stored in
  `diagnostics$solver$cost_scale`.

- rho_query:

  Query-side unbalanced mass relaxation.

- rho_ref:

  Reference-side unbalanced mass relaxation.

- solver:

  Sinkhorn solver variant. `"internal"` (default) and `"auto"` both use
  the primal-domain scaling iteration. `"log_domain"` uses a numerically
  stable log-potential variant that avoids kernel underflow for small
  `epsilon` or high-dimensional codebooks; it is slower per iteration
  but tolerates cost/epsilon ratios that cause `"internal"` to warn.

- min_match_fraction:

  Minimum transported fraction required before a query node label
  transfer is accepted.

- confidence_threshold:

  Minimum top-label probability required before a query node label
  transfer is accepted.

- correction_min_mass:

  Minimum transported node mass required before a correction shift is
  applied. Corrections also require the node match fraction to pass
  `min_match_fraction`.

- max_iter:

  Maximum internal Sinkhorn iterations.

- tol:

  Internal Sinkhorn convergence tolerance.

- chunk_size:

  Integer. Number of samples to project per chunk when computing nearest
  reference node. Use `Inf` or `NULL` for no chunking (allocates a full
  n_samples x n_nodes matrix). Default `10000L`.

## Value

A `somalign_fit` object.

## Details

The transport plan row sums will not equal `query$node_masses` exactly –
this is by design. Unbalanced optimal transport allows mass destruction,
so some query mass may be absorbed rather than transported. Deviation
grows with lower `rho_query` / `rho_ref` values and higher `epsilon`.
Use `diagnostics$ot$max_row_mass_error` to quantify the deviation in a
given fit; for near-balanced data, increase `rho_query` (e.g.
`rho_query = 10`) to enforce tighter marginal constraints. A warning is
emitted automatically when more than 50% of query mass is destroyed.

The cost matrix is normalised by its median positive entry before the
Sinkhorn kernel is computed. This makes `epsilon` scale- and
dimension-invariant: the same value produces the same degree of
regularisation regardless of the number of features or the spread of
codebook coordinates. The raw normalisation factor is stored as
`diagnostics$solver$cost_scale`.

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
```
