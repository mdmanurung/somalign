# Divide each feature column by its upper quantile

Divides each feature column by its upper quantile so that the
`probs`-quantile maps to 1.0. Scale-only normalisation. Returned matrix
passes directly to
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

An optional pre-processing step that brings each marker's dynamic range
into \\\[0, 1\]\\ before projection to the reference. Unlike
[`somalign_normalize()`](https://mdmanurung.github.io/somalign/reference/somalign_normalize.md),
this operates entirely in raw (unscaled) space and does not use the
reference mean or standard deviation; it only uses the reference to
resolve feature names.

## Usage

``` r
somalign_quantile_normalize(data, reference, probs = 0.999, features = NULL)
```

## Arguments

- data:

  Numeric matrix of query data, same format as the `data` argument to
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

- reference:

  A `somalign_reference` object.

- probs:

  Single numeric in the open interval (0, 1). The quantile used as the
  normalisation denominator. Default is 0.999 (99.9th percentile).

- features:

  Optional character vector of feature names. Defaults to
  `reference$features`.

## Value

A numeric matrix with the same dimensions as `data`, with columns in
`reference$features` order. Each column is divided by its
`probs`-quantile so the bulk of the signal maps into \\\[0, 1\]\\. Pass
this matrix directly to
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

## See also

[`somalign_normalize()`](https://mdmanurung.github.io/somalign/reference/somalign_normalize.md),
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)

## Examples

``` r
set.seed(1)
mat <- matrix(abs(rnorm(40)) * 1000, nrow = 20, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
normed <- somalign_quantile_normalize(mat, ref, probs = 0.999)
```
