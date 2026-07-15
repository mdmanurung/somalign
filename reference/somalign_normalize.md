# Global pre-correction of query data to match the reference distribution

An optional pre-processing step that shifts (or shifts and scales) query
data in reference-scaled coordinate space so that the per-marker query
means align with the reference coordinate origin. Passing the returned
matrix to
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
trains the query SOM on pre-centred data, reducing the global component
of the batch shift that
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
then needs to resolve via optimal transport.

## Usage

``` r
somalign_normalize(
  data,
  reference,
  method = c("mean", "scale"),
  features = NULL
)
```

## Arguments

- data:

  Numeric matrix of query data, same format as the `data` argument to
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

- reference:

  A `somalign_reference` object.

- method:

  Normalisation method. `"mean"` (default) subtracts the per-marker
  query mean in reference-scaled space, removing a uniform location
  shift. `"scale"` additionally divides by the per-marker query standard
  deviation, removing a uniform scale shift as well.

- features:

  Optional character vector of feature names. Defaults to
  `reference$features`.

## Value

A numeric matrix with the same number of rows as `data` and columns in
`reference$features` order, expressed in the original (unscaled) units
of `data`. Pass this matrix directly to
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

## Details

`somalign_normalize()` applies the same per-marker shift (and optionally
rescaling) to every cell, so population-specific batch effects remain
for
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
to resolve.

**When not to use this function.** Mean-normalisation assumes the
apparent per-marker shift reflects instrument drift or reagent-lot
differences affecting all populations uniformly. If the shift reflects
genuine compositional differences between batches (e.g.\\ different
cell-type frequencies), subtracting the global mean distorts the
biology.

## See also

[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md),
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
shifted <- mat + 0.5
corrected <- somalign_normalize(shifted, ref)
```
