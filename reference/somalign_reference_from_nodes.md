# Build a reference object from saved node-level artifacts

Build a reference object from saved node-level artifacts

## Usage

``` r
somalign_reference_from_nodes(
  codebook,
  features,
  center,
  scale,
  node_masses = NULL,
  label_prob = NULL,
  distance_quantiles = NULL,
  global_distance_quantiles = NULL
)
```

## Arguments

- codebook:

  Reference node codebook matrix.

- features:

  Feature names and order.

- center:

  Saved reference feature centers.

- scale:

  Saved reference feature scales.

- node_masses:

  Optional reference node masses.

- label_prob:

  Optional node-by-label probability matrix.

- distance_quantiles:

  Optional node-by-quantile distance matrix.

- global_distance_quantiles:

  Optional global reference distance quantiles.

## Value

A `somalign_reference` object.

## Examples

``` r
if (FALSE) { # \dontrun{
cb <- matrix(c(0.1, 0.2, -0.1, 0.3, 0.4, -0.2, 0.0, 0.1),
             nrow = 4, ncol = 2,
             dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_reference_from_nodes(
  codebook = cb,
  features = c("F1", "F2"),
  center   = c(F1 = 0, F2 = 0),
  scale    = c(F1 = 1, F2 = 1)
)
} # }
```
