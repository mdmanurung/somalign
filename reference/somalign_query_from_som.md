# Build a query object from a pre-trained kohonen SOM

Creates a `somalign_query` by reusing the per-cell node assignments
(`som$unit.classif`) already computed during SOM training, bypassing the
O(N \\\times\\ nodes) per-cell argmax that
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
would otherwise perform.

## Usage

``` r
somalign_query_from_som(
  som,
  data,
  reference,
  codebook = NULL,
  codebook_space = c("reference_scaled", "raw"),
  features = NULL
)
```

## Arguments

- som:

  A trained kohonen SOM (output of
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html) or
  [`kohonen::supersom()`](https://rdrr.io/pkg/kohonen/man/supersom.html))
  with `$unit.classif` populated (`keep.data = TRUE` during training is
  *not* required here, but `unit.classif` must be present).

- data:

  Numeric matrix of raw query cell values (cells \\\times\\ features).
  Must include all features in `reference$features`, and `nrow(data)`
  must equal `length(som$unit.classif)`.

- reference:

  A `somalign_reference` object.

- codebook:

  Optional numeric matrix of SOM codebook vectors in the coordinate
  space given by `codebook_space` (nodes \\\times\\ features). When
  `NULL` (default) the X-layer codebook `som$codes[[1]]` is used. Supply
  an explicitly transformed codebook (e.g.\\ after winsorisation and
  rescaling into reference-scaled space) when the SOM's native codebook
  has been post-processed before alignment.

- codebook_space:

  Coordinate system of the codebook. `"reference_scaled"` (default)
  assumes the codebook is already in the reference-scaled coordinate
  system. `"raw"` applies `reference$center` and `reference$scale`
  before use.

- features:

  Optional character vector of feature names. Defaults to
  `reference$features`.

## Value

A `somalign_query` object.

## Details

`somalign_query_from_som()` reuses `som$unit.classif` directly as the
per-cell query-node assignment, so the O(N \\\times\\ nodes) nearest-
code search is skipped entirely. This is exact when the supplied
`codebook` was used as-is during SOM training. If the codebook has been
post-processed (e.g.\\ winsorised and rescaled into reference space),
the reused assignments are an approximation: cells near node boundaries
may flip under the non-linear transform, but in practice this affects
only a small fraction of cells at the distributional tails and does not
measurably impact OT alignment quality.

Unlike
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md),
`sample_distance` is not recomputed and is set to `NA` in the returned
object. The field is not used by
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

The returned object is identical in structure to a `somalign_query`
produced by
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md),
and is fully compatible with
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## See also

[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md),
[`somalign_reference_from_som()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_som.md),
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
