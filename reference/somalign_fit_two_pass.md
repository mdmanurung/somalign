# Two-pass alignment decomposing correction into global and local components

A two-stage variant of
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
that separates the batch correction into a global shift (estimated at
high regularisation) and a local residual (refined at lower
regularisation). This decomposition is most useful when the batch effect
has both a large uniform component and smaller population-specific
residuals.

## Usage

``` r
somalign_fit_two_pass(
  query,
  reference,
  epsilon_global = 0.3,
  epsilon_local = 0.1,
  rho_query = 1,
  rho_ref = 1,
  solver = c("internal", "log_domain", "auto"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L,
  label_guided = FALSE,
  variance_threshold = 0.9
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon_global:

  Entropic regularisation for pass 1 (global). Higher values give a
  smoother, more diffuse transport plan that captures the mean batch
  shift while averaging out population-specific noise. Default `0.3`;
  should be larger than `epsilon_local`.

- epsilon_local:

  Entropic regularisation for pass 2 (local). Should be smaller than
  `epsilon_global` to refine residual node-level corrections. Default
  `0.1`.

- rho_query:

  Query-side unbalanced mass relaxation (both passes).

- rho_ref:

  Reference-side unbalanced mass relaxation (both passes).

- solver:

  Sinkhorn solver variant. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- min_match_fraction:

  Minimum transported fraction for corrections and label transfer. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- confidence_threshold:

  Minimum label confidence for transfer acceptance. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- correction_min_mass:

  Minimum transported mass for correction. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- max_iter:

  Maximum Sinkhorn iterations per pass.

- tol:

  Sinkhorn convergence tolerance.

- chunk_size:

  Integer. Number of samples per projection chunk. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- label_guided:

  Logical. When `TRUE`, applies a large cost penalty to node pairs with
  discordant dominant labels in both OT passes. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  for details.

- variance_threshold:

  Numeric in (0, 1\]. Cumulative singular-value-squared fraction used to
  select the rank of the batch-subspace *diagnostic* stored in
  `$two_pass$batch_subspace`. Default `0.9`. Has no effect on the
  correction.

## Value

A `somalign_fit` object with an additional `$two_pass` list containing
`global_shift` (per-feature vector), `global_shift_norm` (Euclidean
magnitude), `epsilon_global`, `epsilon_local`, and `batch_subspace` (a
list with `V`, `rank`, `variance_explained` derived from the pass-1
correction field — **descriptive only**, not used for correction; may
conflate batch effects with biology).

## Details

Pass 1 runs OT at `epsilon_global` between the original query codebook
and the reference codebook. The mass-weighted mean node shift across
correction-allowed nodes becomes the global shift `g`. Pass 2 runs OT at
`epsilon_local` between the globally shifted query codebook and the
reference, capturing residual population-specific displacements. The
final per-node correction is the residual plus `g`, so the total
correction for each cell equals its pass-2 barycentric target minus its
original codebook centroid.

Direct projection (`old_som_unit`, `old_som_label`, `final_status`) is
computed from the original unshifted `query$scaled_data` and is
unaffected by the transport, preserving the transport-free primary
result that
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
guarantees.

When the batch shift is predominantly global,
`fit$two_pass$global_shift` approximates the per-feature batch offset.
When the shift is negligible, the two-pass result converges toward a
plain
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
at `epsilon_local`.

## See also

[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md),
[`somalign_normalize()`](https://mdmanurung.github.io/somalign/reference/somalign_normalize.md)

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
fit2 <- somalign_fit_two_pass(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.34); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.09); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
```
