# Soft (probabilistic) label projection for query cells

Projects each query cell onto the reference by a Gaussian kernel over
its k nearest reference SOM nodes, and returns a per-cell probability
distribution over labels (or any node-level grouping). This is the soft
analogue of the hard nearest-node assignment behind `old_som_label` in
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md).

## Usage

``` r
somalign_soft_labels(
  fit,
  node_groups = NULL,
  k = 8L,
  bandwidth = NULL,
  chunk_size = 10000L
)
```

## Arguments

- fit:

  A `somalign_fit` object.

- node_groups:

  Optional node-level grouping to project onto instead of the reference
  labels. Either a length-`n_nodes` vector (one group per reference
  node, e.g. a node-to-metacluster map; converted to indicators) or an
  `n_nodes` by `n_groups` matrix of node-group memberships. When `NULL`
  (default), `fit$reference$label_prob` is used, and the reference must
  carry labels.

- k:

  Integer. Number of nearest reference nodes used for the kernel,
  clamped to the number of reference nodes. Default `8L`.

- bandwidth:

  Positive scalar or `NULL`. Gaussian kernel bandwidth in
  reference-scaled space. `NULL` (default) uses the median
  nearest-neighbour distance of the reference codebook.

- chunk_size:

  Positive integer. Cells are processed in blocks of this size to bound
  peak memory. Default `10000L`.

## Value

A numeric matrix of class `c("somalign_soft_labels", "matrix")`, one row
per query cell and one column per label/group, with rows summing to 1 (a
row is all-zero when a cell's nearest nodes carry no label mass). Row
names are the query sample identifiers; attributes `k` and `bandwidth`
record the settings used.

## Details

Hard projection assigns each cell to a single nearest node and inherits
that node's label, discarding the cell's position within the node's
Voronoi region. At a label boundary this makes assignment a
discontinuous 0/1 decision, so a small batch shift can flip a cell's
label and a boundary cell contributes a full unit to one label with no
hedging. Soft projection instead spreads each cell over its nearest
nodes' labels, which removes that boundary discontinuity and reduces the
quantisation variance in downstream per-sample frequency estimates. It
changes the *frequency estimate*, not the most-likely label:
[`max.col()`](https://rdrr.io/r/base/maxCol.html) of a soft-label matrix
typically matches the hard label.

## Coverage contract

Soft projection applies **no** acceptance or out-of-reference gating.
Unlike the hard path in
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
(which flags cells via `outside_reference_distance` /
`transferred_label_accepted` and can suppress low-confidence transfers),
every cell contributes to the soft distribution according to its
distance to the reference nodes alone. Cells that fall outside the
reference are still projected onto their nearest nodes. If you need to
exclude out-of-reference cells, filter them with
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
before aggregating, or subset the query.

## See also

[`somalign_soft_frequencies()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_frequencies.md),
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
  soft <- somalign_soft_labels(fit)
  head(soft)
}
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 7 query node(s) have match_mass_ratio > 1 (max 1.18); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#>           high low
#> 1 2.378258e-34   1
#> 2 6.656011e-25   1
#> 3 1.619808e-28   1
#> 4 6.683424e-28   1
#> 5 9.064764e-25   1
#> 6 2.103095e-30   1
```
