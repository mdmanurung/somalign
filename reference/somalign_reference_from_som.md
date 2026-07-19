# Build a reference object directly from a trained kohonen SOM

Constructs a `somalign_reference` without reprojecting any cells by
reusing information already stored inside the trained kohonen object:

## Usage

``` r
somalign_reference_from_som(
  som,
  center,
  scale,
  codebook_space,
  labels = c("codebook", "none"),
  quantile_probs = c(0.5, 0.9, 0.95, 0.99),
  distance_chunk_size = 1000000L,
  compute_node_var = TRUE
)
```

## Arguments

- som:

  A trained `kohonen` SOM object (from
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html),
  [`kohonen::supersom()`](https://rdrr.io/pkg/kohonen/man/supersom.html),
  or [`kohonen::xyf()`](https://rdrr.io/pkg/kohonen/man/supersom.html))
  with `keep.data = TRUE` (the kohonen default).

- center:

  Named numeric vector of reference feature centers, one per feature.
  Required.

- scale:

  Named numeric vector of reference feature scales, one per feature
  (must be strictly positive). Required.

- codebook_space:

  Coordinate system of the SOM codebook. `"reference_scaled"` when the
  SOM was trained on data already transformed by `center` and `scale`;
  `"raw"` when the SOM was trained on raw feature values that should be
  transformed into reference-scaled space before use.

- labels:

  `"codebook"` (default) reads the per-node label distribution from the
  Y-layer `codes[[2]]` and enables label transfer. `"none"` skips label
  extraction and disables label transfer regardless of whether a Y-layer
  is present.

- quantile_probs:

  Quantile levels for per-node distance thresholds. Passed to
  [`somalign_reference_from_nodes`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md).

- distance_chunk_size:

  Number of cells to process per chunk when computing X-space
  cell-to-node distances. Reduce if memory is tight; increase for faster
  throughput. Default 1e6.

- compute_node_var:

  Logical; if `TRUE` (default) per-node per-marker variances are
  computed from the embedded SOM training data and stored as
  `reference$node_var`. See
  [`somalign_train_reference`](https://mdmanurung.github.io/somalign/reference/somalign_train_reference.md).

## Value

A `somalign_reference` object, identical in structure to the output of
[`somalign_reference`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
but built without reprojecting any cells.

## Details

- X codebook:

  `codes[[1]]` — the reference node positions in feature space, already
  used by the existing API.

- Node masses:

  `tabulate(som$unit.classif)` — exact counts over *all* training cells,
  zero cost.

- Label probabilities:

  `codes[[2]]` — the supervised Y-layer codebook from an
  `xyf`/`supersom` object. Each row is a per-node distribution over
  cell-type labels. Absent for plain
  [`som()`](https://rdrr.io/pkg/kohonen/man/supersom.html) objects, in
  which case label transfer is disabled.

- Distance quantiles:

  Recomputed in reference-scaled X-space from the embedded
  `som$data[[1]]` using each cell's known `unit.classif` assignment.
  This is O(N \\\times\\ p) with no argmax and no O(N \\\times\\ nodes)
  memory peak.

**Note on partition semantics.** On an `xyf`/`supersom` trained with
equal layer weights, `unit.classif` reflects a joint X+Y
(label-weighted) assignment, not a pure-X nearest-node assignment. Node
masses therefore match the SOM's own supervised partition rather than
somalign's X-only geometry. Distance quantiles are still computed in
X-only space so the outside-reference threshold is on the same scale as
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)'s
query distances.

## See also

[`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md),
[`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md),
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md),
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)

## Examples

``` r
set.seed(1)
n <- 120
X <- matrix(rnorm(n * 2), nrow = n, ncol = 2,
            dimnames = list(NULL, c("F1", "F2")))
Y <- cbind(A = rep(c(1, 0), each = n / 2),
           B = rep(c(0, 1), each = n / 2))
som_obj <- kohonen::supersom(list(X, Y),
                             grid = kohonen::somgrid(3, 3, "hexagonal"),
                             rlen = 5, keep.data = TRUE)
center <- colMeans(X)
scale  <- apply(X, 2, sd)
X_scaled <- scale(X, center = center, scale = scale)
som_scaled <- kohonen::supersom(list(X_scaled, Y),
                                grid = kohonen::somgrid(3, 3, "hexagonal"),
                                rlen = 5, keep.data = TRUE)
ref <- somalign_reference_from_som(som_scaled,
                                   center = center, scale = scale,
                                   codebook_space = "reference_scaled")
```
