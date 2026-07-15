# Prepare query data and attach or train a query SOM

Query data are always scaled with the saved reference center and scale.

## Usage

``` r
somalign_query(
  data,
  reference,
  som_query = NULL,
  codebook_space = c("reference_scaled", "raw"),
  grid = NULL,
  rlen = 100,
  alpha = c(0.05, 0.01),
  features = NULL,
  ...
)
```

## Arguments

- data:

  Numeric query data.

- reference:

  A `somalign_reference` object.

- som_query:

  Optional query SOM or SOM-like object with a codebook.

- codebook_space:

  Coordinate system of the `som_query` codebook. Only used when
  `som_query` is supplied. `"reference_scaled"` (default) assumes the
  codebook was already trained on query data scaled with
  `reference$center` and `reference$scale`; `"raw"` re-scales the
  codebook into reference-scaled space before use.

- grid:

  Optional
  [`kohonen::somgrid()`](https://rdrr.io/pkg/kohonen/man/unit.distances.html)
  object when `som_query` is omitted.

- rlen:

  Number of SOM training iterations passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- alpha:

  Learning-rate schedule passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- features:

  Optional feature names. Defaults to the reference feature order.

- ...:

  Additional arguments passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

## Value

A `somalign_query` object.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
```
