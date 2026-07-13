# Build a reference object from an existing SOM and old data

Build a reference object from an existing SOM and old data

## Usage

``` r
somalign_reference(
  som_ref,
  data,
  labels = NULL,
  features = NULL,
  center = NULL,
  scale = NULL,
  codebook_space = NULL,
  quantile_probs = c(0.5, 0.9, 0.95, 0.99)
)
```

## Arguments

- som_ref:

  A `kohonen` SOM object, or a SOM-like object containing a codebook
  matrix.

- data:

  Numeric old/reference data used to compute scaling, masses, labels,
  and distance thresholds.

- labels:

  Optional labels, one per row of `data`.

- features:

  Optional feature names to use. Defaults to all columns.

- center:

  Optional saved feature centers. Computed from `data` when omitted.

- scale:

  Optional saved feature scales. Computed from `data` when omitted.

- codebook_space:

  Coordinate system of the existing `som_ref` codebook. Use
  `"reference_scaled"` when the SOM was trained on `data` transformed
  with `center` and `scale`; use `"raw"` when the SOM was trained on raw
  feature values and should be transformed into reference-scaled space.

- quantile_probs:

  Distance quantiles used for outside-reference flags.

## Value

A `somalign_reference` object.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
g <- kohonen::somgrid(2, 2, "hexagonal")
som_obj <- kohonen::som(scale(mat), grid = g, rlen = 5)
ref <- somalign_reference(som_obj, mat, labels = rep(c("A", "B"), each = 5),
                          codebook_space = "reference_scaled")
```
