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
  range: 0.5–3.

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

## Value

A `somalign_anchored_fit` object (also inherits `somalign_fit`).

## Details

**Cost modification.** Let \\C\\ be the M×K codebook distance matrix
normalised by its median positive entry (as in
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)).
The anchor pairs are projected onto the query codebook (old batch) and
reference codebook (new batch), yielding a count matrix \\A\\ where
\\A\_{kl}\\ is the number of anchor pairs whose old measurement maps to
query node \\k\\ and new measurement maps to reference node \\l\\. (The
query SOM was trained on new-batch data, so projecting the old-batch
anchor onto it identifies which query node the anchor occupied before
the batch shift; projecting the new-batch anchor onto the reference SOM
identifies the corresponding reference node after the shift.) The
modified cost is \$\$\tilde{C}\_{kl} = \max\\\bigl(C\_{kl} -
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
pair reduces cost by at most one median-distance unit.

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

- `nodes_covered`:

  Number of query nodes with ≥ 1 anchor pair.

- `coverage_fraction`:

  `nodes_covered / nrow(query$codebook)`.

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
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.63); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
fit$anchors
#> $n_anchors
#> [1] 10
#> 
#> $rho_anchor
#> [1] 1
#> 
#> $nodes_covered
#> [1] 2
#> 
#> $coverage_fraction
#> [1] 0.5
#> 
```
