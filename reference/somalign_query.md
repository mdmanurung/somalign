# Prepare query data and attach or train a query SOM

Query data are always scaled with the saved reference center and scale.

## Usage

``` r
somalign_query(
  data,
  reference,
  som_query = NULL,
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

  Optional query SOM or SOM-like object with a codebook. The codebook
  must be in the reference-scaled feature space, i.e. trained on query
  data transformed with `reference$center` and `reference$scale`.

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
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
} # }
```
