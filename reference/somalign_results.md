# Return per-sample somalign results

Direct reference projection columns are canonical. Corrected projection
columns are auxiliary and should be used for visualisation, annotation,
and triage rather than feature-level differential testing.

## Usage

``` r
somalign_results(
  fit,
  data = NULL,
  outside_pvalue_threshold = NULL,
  include_correction = TRUE
)
```

## Arguments

- fit:

  A `somalign_fit` object.

- data:

  Optional data frame to append after result columns.

- outside_pvalue_threshold:

  Optional numeric in \[0, 1\]. When supplied, adds a boolean
  `outside_reference_pvalue_flag` column, `TRUE` when
  `outside_reference_pvalue < outside_pvalue_threshold`. `NULL`
  (default) omits the column.

- include_correction:

  Logical. When `FALSE`, drops the auxiliary correction columns
  (`corrected_som_unit`, `corrected_som_distance`,
  `corrected_som_distance_threshold`,
  `corrected_outside_reference_distance`, `correction_norm`) for a
  label-transfer-focused result. The corrected coordinates are a
  diagnostic, not a batch-corrected expression product; see the
  package's topology audit for why they can over-merge populations.
  Default `TRUE` (unchanged behaviour).

## Value

A data frame with direct and corrected projection columns, plus
`transferred_label_second`, `transferred_label_second_confidence`, and
`transferred_label_margin` for triaging low-margin label transfers, and
`outside_reference_surprisal`, `outside_reference_pvalue`, and
`outside_reference_top_marker`: a calibrated chi-squared alternative to
`outside_reference_distance` that weights per-marker deviations from a
query cell's assigned reference node by that node's own per-marker
variance (`reference$node_var`). `NA` when the reference was built
without `node_var` (e.g. `compute_node_var = FALSE`, or a reference
built before this feature). The chi-squared calibration assumes a
diagonal-Gaussian node model; it is anti-conservative for heavy-tailed
(e.g. lognormal) marker distributions, but remains a useful *relative*
ranking of anomalous cells regardless.

## Details

`transferred_label` is the query node's top label transfer choice, an
argmax over `correspondence %*% reference$label_prob`. When two labels
receive close probability mass, this argmax can be brittle: a cell can
land on a node whose true dominant population is clean, yet still get
the runner-up label because it narrowly lost the top spot.
`transferred_label_second` and `transferred_label_margin` (top
confidence minus second confidence, from the same query node) expose
this directly, so close calls can be triaged downstream instead of
silently trusted at face value.

## Cross-batch composition and abundance

For comparing cell-type composition or abundance across batches, use the
**direct** projection columns (`old_som_unit`, `old_som_label`), not the
corrected ones. The barycentric correction can over-merge distinct
populations (see
[`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md))
and does not improve, and can worsen, the cross-batch reproducibility of
composition. Quantify abundance **compositionally** rather than on raw
frequencies: apply a centred log-ratio (CLR) transform to per-sample
cluster counts (for example with the `crumblr` package, which also
supplies count-based precision weights for downstream
differential-abundance models). Raw proportions are dominated by the
largest clusters and understate the reproducibility of rarer ones; the
CLR profile is the reproducible quantity. Absolute per-cluster
proportions and node-level frequencies remain approximate and should not
be treated as precise.

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
fit <- somalign_fit(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_results(fit)
#>    sample_id query_som_unit old_som_unit old_som_distance
#> 1          1              3            1       0.83244247
#> 2          2              3            3       0.46039880
#> 3          3              2            2       0.46022700
#> 4          4              4            4       0.14787688
#> 5          5              3            1       0.75830799
#> 6          6              2            2       0.07897186
#> 7          7              1            3       0.26164018
#> 8          8              3            3       0.69506366
#> 9          9              3            3       0.53031798
#> 10        10              2            1       0.23777951
#>    old_som_distance_threshold outside_reference_distance
#> 1                   0.8250290                       TRUE
#> 2                   0.6703518                      FALSE
#> 3                   0.4411642                       TRUE
#> 4                   0.1478769                      FALSE
#> 5                   0.8250290                      FALSE
#> 6                   0.4411642                      FALSE
#> 7                   0.6703518                      FALSE
#> 8                   0.6703518                       TRUE
#> 9                   0.6703518                      FALSE
#> 10                  0.8250290                      FALSE
#>    outside_reference_surprisal outside_reference_pvalue
#> 1                   2.96492628                0.2270777
#> 2                   2.30068114                0.3165290
#> 3                   2.91411201                0.2329210
#> 4                   0.02186757                0.9891258
#> 5                   1.71095932                0.4250792
#> 6                   0.08580386                0.9580053
#> 7                   0.42477726                0.8086504
#> 8                   3.26341345                0.1955955
#> 9                   1.69971026                0.4274769
#> 10                  0.26978329                0.8738106
#>    outside_reference_top_marker      final_status old_som_label
#> 1                            F2 outside_reference          <NA>
#> 2                            F1  inside_reference          <NA>
#> 3                            F1 outside_reference          <NA>
#> 4                            F1  inside_reference          <NA>
#> 5                            F1  inside_reference          <NA>
#> 6                            F2  inside_reference          <NA>
#> 7                            F2  inside_reference          <NA>
#> 8                            F2 outside_reference          <NA>
#> 9                            F2  inside_reference          <NA>
#> 10                           F2  inside_reference          <NA>
#>    old_som_label_confidence corrected_som_unit corrected_som_distance
#> 1                        NA                  1              0.7582309
#> 2                        NA                  3              0.6913405
#> 3                        NA                  2              0.3598947
#> 4                        NA                  4              0.0000000
#> 5                        NA                  1              0.5213334
#> 6                        NA                  2              0.1918948
#> 7                        NA                  3              0.1097766
#> 8                        NA                  3              0.2575890
#> 9                        NA                  3              0.1875863
#> 10                       NA                  1              0.1931799
#>    corrected_som_distance_threshold corrected_outside_reference_distance
#> 1                         0.8250290                                FALSE
#> 2                         0.6703518                                 TRUE
#> 3                         0.4411642                                FALSE
#> 4                         0.1478769                                FALSE
#> 5                         0.8250290                                FALSE
#> 6                         0.4411642                                FALSE
#> 7                         0.6703518                                FALSE
#> 8                         0.6703518                                FALSE
#> 9                         0.6703518                                FALSE
#> 10                        0.8250290                                FALSE
#>    correction_norm transferred_label transferred_label_confidence
#> 1        0.4395087              <NA>                           NA
#> 2        0.4395087              <NA>                           NA
#> 3        0.1185567              <NA>                           NA
#> 4        0.1478769              <NA>                           NA
#> 5        0.4395087              <NA>                           NA
#> 6        0.1185567              <NA>                           NA
#> 7        0.3455139              <NA>                           NA
#> 8        0.4395087              <NA>                           NA
#> 9        0.4395087              <NA>                           NA
#> 10       0.1185567              <NA>                           NA
#>    transferred_label_accepted transferred_label_second
#> 1                       FALSE                     <NA>
#> 2                       FALSE                     <NA>
#> 3                       FALSE                     <NA>
#> 4                       FALSE                     <NA>
#> 5                       FALSE                     <NA>
#> 6                       FALSE                     <NA>
#> 7                       FALSE                     <NA>
#> 8                       FALSE                     <NA>
#> 9                       FALSE                     <NA>
#> 10                      FALSE                     <NA>
#>    transferred_label_second_confidence transferred_label_margin
#> 1                                   NA                       NA
#> 2                                   NA                       NA
#> 3                                   NA                       NA
#> 4                                   NA                       NA
#> 5                                   NA                       NA
#> 6                                   NA                       NA
#> 7                                   NA                       NA
#> 8                                   NA                       NA
#> 9                                   NA                       NA
#> 10                                  NA                       NA
```
