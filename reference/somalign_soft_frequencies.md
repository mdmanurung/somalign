# Per-group soft label frequencies for query cells

Aggregates the per-cell soft label distributions from
[`somalign_soft_labels()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_labels.md)
by a grouping (typically a biological sample), giving a group-by-label
matrix of soft frequencies. Soft aggregation reduces the
boundary/quantisation variance of hard per-sample cluster proportions,
which improves the reproducibility of cluster-abundance profiles across
batches (for example the centred-log-ratio abundance comparison in the
label-transfer vignette).

## Usage

``` r
somalign_soft_frequencies(
  fit,
  group,
  node_groups = NULL,
  k = 8L,
  bandwidth = NULL,
  normalize = TRUE,
  chunk_size = 10000L
)
```

## Arguments

- fit:

  A `somalign_fit` object.

- group:

  Vector of length equal to the number of query cells, giving each
  cell's group (e.g. `sample_id` or `fcs_filename`).

- node_groups:

  Optional node-level grouping passed to
  [`somalign_soft_labels()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_labels.md)
  (e.g. a node-to-metacluster map). Default `NULL` uses the reference
  labels.

- k, bandwidth, chunk_size:

  Passed to
  [`somalign_soft_labels()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_labels.md).

- normalize:

  Logical. When `TRUE` (default) each group's row is divided by its
  total so rows are frequencies summing to 1; when `FALSE` the raw
  summed soft memberships (soft counts) are returned, suitable for
  count-based differential-abundance models. Note that a group's soft
  counts sum to the number of that group's cells whose neighbours carry
  label mass, not necessarily its total cell count: cells all of whose k
  nearest nodes are unlabelled contribute a zero row. With a fully
  labelled reference the two coincide.

## Value

A numeric matrix of class `c("somalign_soft_frequencies", "matrix")`,
one row per group and one column per label/group, with attributes `k`,
`bandwidth`, and `normalized`.

## See also

[`somalign_soft_labels()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_labels.md),
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)

## Examples

``` r
if (requireNamespace("kohonen", quietly = TRUE)) {
  set.seed(1)
  x <- rbind(matrix(rnorm(90 * 3, -3, 0.5), ncol = 3),
             matrix(rnorm(90 * 3,  3, 0.5), ncol = 3))
  colnames(x) <- paste0("m", seq_len(3))
  lab <- rep(c("low", "high"), each = 90)
  grid <- kohonen::somgrid(3, 3, "hexagonal")
  ref <- somalign_train_reference(x, labels = lab, grid = grid, rlen = 10)
  qry <- somalign_query(x, ref, grid = grid, rlen = 10)
  fit <- somalign_fit(qry, ref)
  sample_id <- rep(c("s1", "s2", "s3"), length.out = nrow(x))
  somalign_soft_frequencies(fit, sample_id)
}
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 7 query node(s) have match_mass_ratio > 1 (max 1.18); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> <somalign_soft_frequencies> [3 groups x 2 labels]  frequencies  k = 8  bandwidth = 0.2356
```
