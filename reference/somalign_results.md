# Return per-sample somalign results

Direct reference projection columns are canonical. Corrected projection
columns are auxiliary and should be used for visualisation, annotation,
and triage rather than feature-level differential testing.

## Usage

``` r
somalign_results(fit, data = NULL)
```

## Arguments

- fit:

  A `somalign_fit` object.

- data:

  Optional data frame to append after result columns.

## Value

A data frame with direct and corrected projection columns.

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
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.18); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
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
#>    old_som_distance_threshold outside_reference_distance      final_status
#> 1                   0.8250290                       TRUE outside_reference
#> 2                   0.6703518                      FALSE  inside_reference
#> 3                   0.4411642                       TRUE outside_reference
#> 4                   0.1478769                      FALSE  inside_reference
#> 5                   0.8250290                      FALSE  inside_reference
#> 6                   0.4411642                      FALSE  inside_reference
#> 7                   0.6703518                      FALSE  inside_reference
#> 8                   0.6703518                       TRUE outside_reference
#> 9                   0.6703518                      FALSE  inside_reference
#> 10                  0.8250290                      FALSE  inside_reference
#>    old_som_label old_som_label_confidence corrected_som_unit
#> 1           <NA>                       NA                  1
#> 2           <NA>                       NA                  3
#> 3           <NA>                       NA                  2
#> 4           <NA>                       NA                  4
#> 5           <NA>                       NA                  1
#> 6           <NA>                       NA                  2
#> 7           <NA>                       NA                  3
#> 8           <NA>                       NA                  3
#> 9           <NA>                       NA                  3
#> 10          <NA>                       NA                  1
#>    corrected_som_distance corrected_som_distance_threshold
#> 1            7.792775e-01                        0.8250290
#> 2            7.047406e-01                        0.6703518
#> 3            4.054800e-01                        0.4411642
#> 4            5.641698e-06                        0.1478769
#> 5            5.009436e-01                        0.8250290
#> 6            1.656497e-01                        0.4411642
#> 7            1.097853e-01                        0.6703518
#> 8            2.627448e-01                        0.6703518
#> 9            2.081463e-01                        0.6703518
#> 10           2.406763e-01                        0.8250290
#>    corrected_outside_reference_distance correction_norm transferred_label
#> 1                                 FALSE       0.4390426              <NA>
#> 2                                  TRUE       0.4390426              <NA>
#> 3                                 FALSE       0.1056234              <NA>
#> 4                                 FALSE       0.1478819              <NA>
#> 5                                 FALSE       0.4390426              <NA>
#> 6                                 FALSE       0.1056234              <NA>
#> 7                                 FALSE       0.3455129              <NA>
#> 8                                 FALSE       0.4390426              <NA>
#> 9                                 FALSE       0.4390426              <NA>
#> 10                                FALSE       0.1056234              <NA>
#>    transferred_label_confidence transferred_label_accepted
#> 1                            NA                      FALSE
#> 2                            NA                      FALSE
#> 3                            NA                      FALSE
#> 4                            NA                      FALSE
#> 5                            NA                      FALSE
#> 6                            NA                      FALSE
#> 7                            NA                      FALSE
#> 8                            NA                      FALSE
#> 9                            NA                      FALSE
#> 10                           NA                      FALSE
```
