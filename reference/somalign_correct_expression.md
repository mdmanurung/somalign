# Batch-correct query marker expression for downstream analysis

Returns a cell-level (cells by markers) batch-corrected marker
expression matrix for the query cells, intended for downstream
visualisation and differential expression. The correction is restricted
to the anchor-estimated batch subspace and smoothed across each cell's
nearest self-organising map (SOM) nodes, so variation orthogonal to the
batch direction is preserved.

## Usage

``` r
somalign_correct_expression(
  fit,
  units = c("raw", "scaled"),
  smooth = TRUE,
  k = 8L,
  bandwidth = NULL,
  confidence_gate = TRUE,
  chunk_size = 10000L
)
```

## Arguments

- fit:

  A `somalign_fit` from
  [`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
  with `correction = "subspace"` or `"both"`, or from
  [`somalign_fit_two_pass()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_two_pass.md).
  A plain
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  or a `correction = "cost_bonus"` anchored fit carries no batch
  subspace and is rejected.

- units:

  One of `"raw"` (original expression units, the default) or `"scaled"`
  (reference-scaled units).

- smooth:

  Logical. When `TRUE` (default), smooths the correction across the k
  nearest SOM nodes with a Gaussian kernel. When `FALSE`, each cell
  takes its nearest node's shift directly (piecewise constant); this
  reproduces the correction `somalign` uses internally and is provided
  as a diagnostic baseline, not recommended for downstream analysis.

- k:

  Integer. Number of nearest SOM nodes used for smoothing, clamped to
  the number of query SOM nodes. Default `8L`.

- bandwidth:

  Positive scalar or `NULL`. Gaussian kernel bandwidth in
  reference-scaled space. `NULL` (default) uses the median
  nearest-neighbour distance of the SOM codebook, which adapts to the
  lattice spacing.

- confidence_gate:

  Logical. When `TRUE` (default), each node's kernel weight is
  multiplied by its transported match fraction, down-weighting nodes the
  transport plan could not align. Nodes whose correction is disallowed
  contribute zero weight in either case. Has no effect when
  `smooth = FALSE`, where each cell simply takes its nearest node's
  shift with no kernel weighting.

- chunk_size:

  Positive integer. Cells are processed in blocks of this size to bound
  peak memory. Default `10000L`.

## Value

A numeric matrix of class
`c("somalign_corrected_expression", "matrix")`, with one row per query
cell and one column per marker. Row names are the query sample
identifiers and column names are the reference features. Attributes
`units`, `bandwidth`, `smooth`, and `k` record the settings used.

## Details

Unlike the per-node correction used internally (which is piecewise
constant across a node and contracts populations toward one another),
this function interpolates a smooth per-cell shift from the shifts of
the k nearest SOM nodes and confines it to the batch subspace. Cells
therefore receive a continuous correction, and structure orthogonal to
the batch subspace is left intact.

## Scope and limitations

This output is an auxiliary correction aid, not the primary product of
`somalign`, which is label transfer (see
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)).
The correction is restricted to the batch subspace and preserves
orthogonal variation, but within that subspace it still *reduces*,
rather than fully removes, the distance between populations; it does not
undo genuine over-merging. For comparing cell-type composition or
abundance across batches, use the direct projection columns from
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
together with a compositional (centred log-ratio) transform, not
corrected expression. Run
[`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md)
before relying on this output to confirm that correction is warranted
for your data.

## Subspace restriction

For an anchored `correction = "subspace"`/`"both"` fit the node shifts
are projected onto the span of the batch subspace \\V\\ at fit time.
Each cell shift is a weighted average of those node shifts, so it lies
in span(\\V\\) too (no post-smoothing re-projection is applied), and
variation orthogonal to \\V\\ is untouched. A two-pass fit does **not**
carry an anchor-estimated subspace: its stored shifts are the full
documented two-pass correction (population-specific residual plus global
shift), and `somalign_correct_expression()` applies them in full without
subspace confinement. The `$two_pass$batch_subspace` diagnostic is not
used here.

## Future extension

A contraction-free variant would estimate the batch-shift field directly
from anchor displacements by kernel regression, rather than from the
barycentric node shifts. That path needs the scaled anchor positions,
which a fit does not store (only the anchor displacement matrix is
retained), so it is not available here.

## See also

[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md),
[`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md),
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md),
[`somalign_fit_two_pass()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_two_pass.md)

## Examples

``` r
if (requireNamespace("kohonen", quietly = TRUE)) {
  set.seed(1)
  ref_x <- matrix(rnorm(60 * 3, 0, 0.5), ncol = 3,
                  dimnames = list(NULL, paste0("m", seq_len(3))))
  grid <- kohonen::somgrid(2, 2, "hexagonal")
  ref <- somalign_train_reference(ref_x, grid = grid, rlen = 10)
  shift <- matrix(c(2, 0, 0), nrow(ref_x), 3, byrow = TRUE)
  qry <- somalign_query(ref_x + shift, ref, grid = grid, rlen = 10)
  anc <- ref_x[1:20, ]
  fit <- somalign_fit_anchored(qry, ref, anchor_old = anc,
                               anchor_new = anc + shift[1:20, ],
                               correction = "subspace")
  expr <- somalign_correct_expression(fit)
  dim(expr)
}
#> Warning: 98.3% of query samples project outside reference distance thresholds. This may indicate a distributional mismatch or a coordinate-space misconfiguration.
#> [1] 60  3
```
