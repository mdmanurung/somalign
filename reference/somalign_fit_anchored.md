# Align a query SOM to a reference SOM using anchor sample pairs

A variant of
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
for the case where a set of samples has been measured in **both** the
old batch (reference space) and the new batch (query space). These
*anchor pairs* are used to build a per-node-pair correspondence count
matrix, which is subtracted from the normalized OT cost before the
Sinkhorn solve. This makes transport along anchor-supported routes
cheaper, biasing the OT plan toward pairings that are consistent with
the observed per-sample batch displacement — while still solving a valid
optimal transport problem over the full codebook.

## Usage

``` r
somalign_fit_anchored(
  query,
  reference,
  anchor_old,
  anchor_new,
  rho_anchor = 1,
  epsilon = 0.1,
  rho_query = 1,
  rho_ref = 1,
  solver = c("internal", "log_domain", "auto", "annealing"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L,
  correction = c("cost_bonus", "subspace", "both"),
  variance_threshold = 0.9,
  anneal_start = 10,
  anneal_stages = 10L,
  anneal_factor = NULL,
  feature_weights = NULL,
  laplacian_lambda = 0
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- anchor_old:

  Numeric matrix (n_anchors × p). Old-batch measurements of the anchor
  samples. Must be **raw (un-normalized) values in the same units and
  preprocessing pipeline as the data used to train `reference`**. Do not
  pre-center or pre-scale; this function applies `reference$center` and
  `reference$scale` internally. Also accepts a data frame of numeric
  columns.

- anchor_new:

  Numeric matrix (n_anchors × p). New-batch measurements of the **same**
  anchor samples. Must be **raw (un-normalized) values in the same units
  and preprocessing pipeline as `anchor_old`**. Do not pre-center or
  pre-scale; this function applies `reference$center` and
  `reference$scale` internally. Rows of `anchor_old` and `anchor_new`
  must correspond to the same biological units. Also accepts a data
  frame of numeric columns.

- rho_anchor:

  Non-negative scalar. Controls how strongly anchor pairs bias the OT
  cost. At `rho_anchor = 0` the result equals
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  Larger values reduce the effective cost for anchor-supported node
  pairs, concentrating the transport plan on those routes. Typical
  range: 0.5–3. Has no effect when `correction = "subspace"`.

- epsilon:

  Entropic regularisation strength (see
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)).

- rho_query:

  Query-side unbalanced mass relaxation.

- rho_ref:

  Reference-side unbalanced mass relaxation.

- solver:

  Sinkhorn solver variant. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- min_match_fraction:

  Minimum transported fraction for label transfer.

- confidence_threshold:

  Minimum top-label probability for label transfer.

- correction_min_mass:

  Minimum transported mass for a node correction.

- max_iter:

  Maximum Sinkhorn iterations.

- tol:

  Sinkhorn convergence tolerance.

- chunk_size:

  Integer. Samples projected per chunk. Default `10000L`.

- correction:

  Character. Correction strategy — one of `"cost_bonus"` (default),
  `"subspace"`, or `"both"`. See Details.

- variance_threshold:

  Numeric in (0, 1\]. Cumulative singular-value-squared fraction for
  selecting the rank of the batch subspace. Default `0.9` (CellANOVA
  convention). Only used when `correction` is `"subspace"` or `"both"`.

- anneal_start, anneal_stages, anneal_factor:

  Annealing-schedule tuning parameters, used only when
  `solver = "annealing"`. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- feature_weights:

  Either `NULL` (default, squared-Euclidean cost), a named non-negative
  numeric vector of explicit per-feature weights (see
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)),
  or the string `"anchor"` – auto-estimates weights from the anchor
  displacement matrix `D` via \\w_f = 1 / (\mathrm{var}(D\_{\cdot f}) +
  \delta)\\, mean-normalised. Markers that vary most across the batch
  (large `var(D[, f])`, i.e. batch-driven) get low weight and are cheap
  to transport; markers stable across the batch get high weight and are
  expensive to transport, preserving biology. The resolved vector is
  stored in `fit$anchors$feature_weights` and
  `fit$diagnostics$cost_metric$feature_weights`. Composes independently
  with `correction`: the weights reshape the cost geometry, while
  `rho_anchor`/`correction` bias routing – both act on the same
  underlying transport problem without conflict.

- laplacian_lambda:

  Non-negative scalar. Graph-Laplacian smoothing of the node-shift
  field; see
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  When `correction` is `"subspace"` or `"both"`, smoothing is applied
  *before* the subspace projection (smooth in full marker space, then
  restrict to the batch subspace `V`) so the Laplacian neighbor
  structure is respected. Default `0` (no smoothing).

## Value

A `somalign_anchored_fit` object (also inherits `somalign_fit`).

## Details

**Correction modes.** Three strategies are available via the
`correction` argument.

- `"cost_bonus"` (default, current behaviour): the anchor count matrix
  biases the OT cost so anchor-supported node pairs are cheaper; the
  resulting node shifts are applied to the full feature space.

- `"subspace"`: a batch subspace \\V\_{\text{batch}}\\ is estimated by
  SVD of the anchor displacement matrix \\D = X\_{\text{old}} -
  X\_{\text{new}}\\ (n_anchors × p). Because each row of \\D\\ is a
  *same-biological-unit* before–after measurement, the dominant singular
  vectors isolate the true batch direction. Node shifts from a *plain*
  OT solve (no cost bonus) are then projected onto
  \\V\_{\text{batch}}\\: only the batch-direction component is applied.
  Biological variation orthogonal to \\V\_{\text{batch}}\\ is preserved.
  A synthetic validation shows the orthogonal component survives at
  ~99.7% (1.496 vs ideal 1.500) while `"cost_bonus"` erases it. The rank
  \\r\\ is the smallest index where the cumulative squared singular
  values reach `variance_threshold` (default 0.9). \\D\\ is **not
  centred** — the mean batch direction is the dominant structure we want
  to capture.

- `"both"`: applies the cost bonus to the OT solve *and* restricts the
  resulting shifts to \\V\_{\text{batch}}\\.

`"subspace"` and `"both"` expose `fit$anchors$batch_subspace` (a list
with `V`, `rank`, `variance_explained`). `"cost_bonus"` sets this to
`NULL`.

**Topology preservation and epsilon.** Empirically, the primary driver
of topology/structure damage from batch correction is `epsilon`, not
`rho_anchor`. Higher epsilon blurs the transport plan across a wider
neighbourhood, causing the corrected codebook to collapse biologically
distinct populations (H0 component merging). Subspace-restricted modes
(`"subspace"` or `"both"`) substantially reduce merging at any given
epsilon because shifts are confined to the batch-variation subspace,
leaving orthogonal biological variation intact. As a result, choosing
epsilon involves a genuine trade-off: higher epsilon is more numerically
stable for the Sinkhorn solver, but lower epsilon preserves more
topology. Before committing to an epsilon, run
`somalign_epsilon_sweep(..., topology = TRUE)` alongside
[`somalign_select_epsilon()`](https://mdmanurung.github.io/somalign/reference/somalign_select_epsilon.md)
and inspect both the phase-transition criterion and the
`biggest_merge_mass_frac` column – the two criteria can disagree,
especially at small epsilon near numerical instability.

**Cost modification.** Let \\C\\ be the M×K codebook distance matrix
normalised by its median positive entry (as in
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)).
Each anchor pair is projected onto both codebooks to build a count
matrix \\A\\ where \\A\_{kl}\\ is the number of anchor pairs mapping to
query node \\k\\ and reference node \\l\\. The query SOM was trained on
new-batch data, so the *new-batch* anchor measurement is projected onto
the query codebook to identify query node \\k\\; the reference SOM was
trained on old-batch data, so the *old-batch* anchor measurement is
projected onto the reference codebook to identify reference node \\l\\.
The modified cost is \$\$\tilde{C}\_{kl} = \max\\\bigl(C\_{kl} -
\rho\_{\mathrm{anchor}} \cdot A\_{kl} / n\_{\mathrm{anchors}},\\
0\bigr).\$\$ Pairs with many anchor observations get cost reduced toward
zero (free transport), while uncovered pairs retain their original cost.
Non-negativity is enforced by the \\\max(\cdot, 0)\\ clamp.

**Clamp behaviour at large `rho_anchor`.** When the anchor bonus exceeds
\\C\_{kl}\\, the effective cost is clamped to zero. All such pairs then
have identical effective cost and the transport mass among them is
determined by entropic regularisation alone rather than by relative
anchor counts. The clamp is required to keep costs non-negative; at very
large `rho_anchor` the plan for anchor-covered pairs becomes more
entropic, not more concentrated. A practical upper bound is
`rho_anchor * max(A) / n_anchors <= 1`, i.e., even the most-supported
pair reduces cost by at most one median squared-distance unit.

**Fallback for uncovered nodes.** Query nodes with no anchor samples
retain their original pairwise costs, so the transport plan for those
nodes is determined entirely by the OT objective — the same as
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
Inspect `$anchors$coverage_fraction` to see what fraction of query nodes
had at least one anchor pair.

**Return value.** The object has class
`c("somalign_anchored_fit", "somalign_fit")`, so all downstream
functions that accept a `somalign_fit` object
([`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md),
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md))
work unchanged. An additional `$anchors` list element is attached:

- `n_anchors`:

  Number of anchor pairs supplied.

- `rho_anchor`:

  The value of `rho_anchor` used.

- `correction`:

  The correction mode: `"cost_bonus"`, `"subspace"`, or `"both"`.

- `nodes_covered`:

  Number of query nodes with ≥ 1 anchor pair.

- `coverage_fraction`:

  `nodes_covered / nrow(query$codebook)`.

- `batch_subspace`:

  For `"subspace"` and `"both"` modes: a list with `V` (p × rank
  matrix), `rank` (integer), and `variance_explained` (cumulative
  variance at the selected rank). `NULL` for `"cost_bonus"`.

- `displacements`:

  The scaled anchor displacement matrix \\D = X\_{\text{old,scaled}} -
  X\_{\text{new,scaled}}\\ (n_anchors × p), always stored regardless of
  `correction` mode. Used by
  [`somalign_subspace_sensitivity()`](https://mdmanurung.github.io/somalign/reference/somalign_subspace_sensitivity.md)
  and
  [`somalign_exclusion_test()`](https://mdmanurung.github.io/somalign/reference/somalign_exclusion_test.md).

## Note

At small `epsilon` with high anchor coverage the anchor bonus zeros out
many entries of the normalised cost matrix, which sharpens the Sinkhorn
kernel and can drive the remaining entries toward numerical underflow.
If the solver warns about kernel underflow, pass `solver = "log_domain"`
or `solver = "annealing"`, both of which work in log-potential space and
avoid the issue.

## See also

[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
for the unanchored variant.

## Examples

``` r
set.seed(1)
p   <- 3L
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
# Use 10 samples as anchors measured in both batches
anc_idx <- 1:10
fit <- somalign_fit_anchored(qry, ref,
                              anchor_old = mat[anc_idx, , drop = FALSE],
                              anchor_new = shifted[anc_idx, , drop = FALSE],
                              rho_anchor = 1)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.11); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
fit$anchors
#> $n_anchors
#> [1] 10
#> 
#> $rho_anchor
#> [1] 1
#> 
#> $correction
#> [1] "cost_bonus"
#> 
#> $nodes_covered
#> [1] 2
#> 
#> $coverage_fraction
#> [1] 0.5
#> 
#> $batch_subspace
#> NULL
#> 
#> $variance_threshold
#> [1] 0.9
#> 
#> $displacements
#>               F1        F2         F3
#>  [1,] -0.2266608 -0.220545 -0.2307808
#>  [2,] -0.2266608 -0.220545 -0.2307808
#>  [3,] -0.2266608 -0.220545 -0.2307808
#>  [4,] -0.2266608 -0.220545 -0.2307808
#>  [5,] -0.2266608 -0.220545 -0.2307808
#>  [6,] -0.2266608 -0.220545 -0.2307808
#>  [7,] -0.2266608 -0.220545 -0.2307808
#>  [8,] -0.2266608 -0.220545 -0.2307808
#>  [9,] -0.2266608 -0.220545 -0.2307808
#> [10,] -0.2266608 -0.220545 -0.2307808
#> 
#> $feature_weights
#> NULL
#> 
```
