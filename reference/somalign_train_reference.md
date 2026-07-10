# Train a reference SOM and build a somalign reference

Train a reference SOM and build a somalign reference

## Usage

``` r
somalign_train_reference(
  data,
  labels = NULL,
  features = NULL,
  grid = NULL,
  rlen = 100,
  alpha = c(0.05, 0.01),
  ...
)
```

## Arguments

- data:

  Numeric matrix or data frame containing old/reference samples.

- labels:

  Optional labels, one per row of `data`.

- features:

  Optional feature names to use. Defaults to all columns.

- grid:

  Optional
  [`kohonen::somgrid()`](https://rdrr.io/pkg/kohonen/man/unit.distances.html)
  object.

- rlen:

  Number of SOM training iterations passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- alpha:

  Learning-rate schedule passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

- ...:

  Additional arguments passed to
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html).

## Value

A `somalign_reference` object.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
labels <- rep(c("A", "B"), each = 5)
ref <- somalign_train_reference(mat, labels = labels,
                                grid = kohonen::somgrid(2, 2, "hexagonal"))
} # }
```
