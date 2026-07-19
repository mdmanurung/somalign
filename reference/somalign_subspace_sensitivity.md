# Subspace sensitivity analysis for anchored batch correction

Bootstraps the anchor displacement matrix `D` to quantify how stable the
estimated batch subspace `V` is, and propagates that uncertainty to
per-node correction confidence intervals and "tipping angles" – the
smallest rotation of `V` that would erase a node's correction.

## Usage

``` r
somalign_subspace_sensitivity(
  fit,
  n_boot = 200L,
  variance_threshold = NULL,
  conf_level = 0.95,
  seed = 1L
)
```

## Arguments

- fit:

  A `somalign_anchored_fit` with `correction = "subspace"` or `"both"`.

- n_boot:

  Positive integer. Number of bootstrap replicates of `D`. Default
  `200L`.

- variance_threshold:

  Numeric in (0, 1\], or `NULL`. Variance threshold for SVD rank
  selection in bootstrap replicates. `NULL` (default) reuses the
  threshold from the original fit (`fit$anchors$variance_threshold`).

- conf_level:

  Numeric in (0, 1). Confidence level for the node-shift confidence
  intervals. Default `0.95`.

- seed:

  Integer or `NULL`. RNG seed for reproducibility (restored on exit;
  does not leak into the caller's session). Default `1L`; `NULL`
  disables seeding.

## Value

A list of class `somalign_subspace_sensitivity`:

- `node_correction_ci`:

  M x 2 matrix (lower, upper) – bootstrap CI on each node's
  corrected-shift norm.

- `node_shift_ci`:

  M x p x 2 array – per-feature bootstrap CI.

- `subspace_angles`:

  n_boot x rank matrix of principal angles (degrees) between the fitted
  `V` and each bootstrap `V_b`.

- `tipping_angle_deg`:

  Length-M vector. For rank-1 `V`, the analytic angle
  \\\arcsin(\|\hat{s}\_i \cdot v\|)\\ between a node's unit raw shift
  and `v`; for rank \> 1, a conservative proxy from the minimum singular
  value of the (1 x rank) projection. `NA` for disallowed or zero-norm
  nodes. Small (\< 10 degrees) signals a fragile correction; large (\>
  45 degrees) indicates robustness.

- `anchor_leverage`:

  Length-n_anchors vector: the maximum principal angle (degrees) between
  `V` and the leave-one-out subspace when that anchor is dropped – a
  Cook's-distance analog for anchor influence.

- `n_boot`, `conf_level`, `subspace_rank`, `n_anchors`,
  `variance_threshold`:

  Metadata.

## Details

Raw pre-projection node shifts are recovered from stored fit components
(`fit$correspondence`, both codebooks), not from `fit$node_shifts`,
which stores *post-projection* shifts (`S \%*\% V \%*\% t(V)`) for
subspace fits – using it directly as the raw shift would double-project
and produce invalid tipping angles.

Bootstrap replicate rank is fixed to the point-estimate rank of the
original fit (not re-selected per replicate), so principal angles
compare subspaces of the same dimension; if a replicate's own rank is
lower, its `V_b` is zero-padded.

## Role

A **diagnostic on the correction path** (how trustworthy the subspace
correction is), not a label-transfer diagnostic. Label transfer does not
use the subspace or `node_shifts`.

## See also

[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md),
[`somalign_exclusion_test()`](https://mdmanurung.github.io/somalign/reference/somalign_exclusion_test.md)

## Examples

``` r
set.seed(1)
p <- 3L
mat <- rbind(
  matrix(rnorm(20 * p, mean = -2), ncol = p),
  matrix(rnorm(20 * p, mean =  2), ncol = p)
)
colnames(mat) <- paste0("F", seq_len(p))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
shifted <- mat + 0.5
qry <- somalign_query(shifted, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
anc_idx <- 1:10
fit <- somalign_fit_anchored(qry, ref,
                              anchor_old = mat[anc_idx, , drop = FALSE],
                              anchor_new = shifted[anc_idx, , drop = FALSE],
                              rho_anchor = 1, correction = "subspace")
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.11); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_subspace_sensitivity(fit, n_boot = 50L)
#> <somalign_subspace_sensitivity>
#>   n_anchors = 10   rank = 1   n_boot = 50   conf_level = 0.95
#>   median principal angle: 0.0 deg
#>   median tipping angle:    79.2 deg (n=4 allowed nodes)
```
