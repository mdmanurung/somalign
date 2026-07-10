# Using trained old and new SOMs

This workflow applies when you have an existing reference SOM and a
separately trained query SOM. Both codebooks must be in the same feature
coordinate system: the reference-scaled space.

## Building a reference from an existing SOM

When building a reference from a saved old SOM, state the codebook
coordinate system explicitly:

``` r

reference <- somalign_reference(
  old_som,
  old_data,
  codebook_space = "reference_scaled"
)
```

Use `codebook_space = "raw"` only if the old SOM codebook is in raw
feature units and should be internally transformed with the reference
center and scale.

## Training the query SOM in reference-scaled space

The query SOM must be trained on new samples transformed with the
**reference** center and scale, not in raw feature units or new-data
z-scores:

``` r

new_scaled_for_som <- sweep(
  sweep(new_data[, reference$features], 2, reference$center, "-"),
  2,
  reference$scale,
  "/"
)
new_som <- kohonen::som(new_scaled_for_som, grid = ...)
```

## Full example

``` r

library(kohonen)
library(somalign)

set.seed(2)
old_data <- rbind(
  matrix(rnorm(30 * 6, mean = -1), ncol = 6),
  matrix(rnorm(30 * 6, mean = 1), ncol = 6)
)
colnames(old_data) <- paste0("marker_", seq_len(ncol(old_data)))
old_labels <- rep(c("old_low", "old_high"), each = 30)

reference <- somalign_train_reference(
  old_data,
  labels = old_labels,
  grid = kohonen::somgrid(3, 2, "hexagonal"),
  rlen = 10
)

new_data <- rbind(
  old_data[1:15, ] + 0.15,
  old_data[46:60, ] + 0.25,
  matrix(rnorm(8 * 6, mean = 4), ncol = 6)
)
colnames(new_data) <- colnames(old_data)

new_scaled_for_som <- sweep(
  sweep(new_data[, reference$features], 2, reference$center, "-"),
  2,
  reference$scale,
  "/"
)

new_som <- kohonen::som(
  new_scaled_for_som,
  grid = kohonen::somgrid(3, 2, "hexagonal"),
  rlen = 10
)

query <- somalign_query(
  new_data,
  reference,
  som_query = new_som
)

fit <- somalign_fit(query, reference, solver = "internal")
results <- somalign_results(fit)

head(results[, c(
  "old_som_unit",
  "old_som_label",
  "outside_reference_distance",
  "final_status",
  "corrected_som_unit",
  "corrected_outside_reference_distance",
  "correction_norm",
  "transferred_label",
  "transferred_label_accepted"
)])
#>   old_som_unit old_som_label outside_reference_distance     final_status
#> 1            1       old_low                      FALSE inside_reference
#> 2            5       old_low                      FALSE inside_reference
#> 3            1       old_low                      FALSE inside_reference
#> 4            1       old_low                      FALSE inside_reference
#> 5            4       old_low                      FALSE inside_reference
#> 6            1       old_low                      FALSE inside_reference
#>   corrected_som_unit corrected_outside_reference_distance correction_norm
#> 1                  1                                FALSE       0.3358327
#> 2                  5                                 TRUE       0.7297865
#> 3                  5                                FALSE       1.2892911
#> 4                  1                                FALSE       0.3358327
#> 5                  4                                FALSE       0.3358327
#> 6                  1                                FALSE       0.3358327
#>   transferred_label transferred_label_accepted
#> 1           old_low                       TRUE
#> 2           old_low                       TRUE
#> 3          old_high                       TRUE
#> 4           old_low                       TRUE
#> 5           old_low                       TRUE
#> 6           old_low                       TRUE
```

## Interpreting the output

**Direct projection columns** (`old_som_unit`, `old_som_distance`,
`outside_reference_distance`, `final_status`, `old_som_label`,
`old_som_label_confidence`) are the primary result. They are based on
nearest-node assignment and are independent of the OT correction.

**OT-corrected columns** (`corrected_som_unit`,
`corrected_som_distance`, `corrected_outside_reference_distance`,
`correction_norm`) are auxiliary.
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
aligns query SOM nodes to reference SOM nodes via OT and computes
node-level correction shifts. For each sample, the shift is applied to
the query SOM position and the corrected position is projected to the
nearest reference node. Use `corrected_som_unit` for visualisation and
triage; keep `old_som_unit` and `final_status` for conservative
classification.

**Label transfer columns** (`transferred_label`,
`transferred_label_confidence`, `transferred_label_accepted`) propagate
reference node labels through the OT correspondence. A transfer is
accepted only when the query node’s match fraction ≥
`min_match_fraction` (default 0.05) and the top-label confidence of the
matched reference node ≥ `confidence_threshold` (default 0.6). Compare
`transferred_label` against `old_som_label` as a cross-check: large
systematic differences indicate a poor OT alignment or a high-novelty
query.

## Checking OT quality

Inspect `diagnostics$ot$match_fraction` before trusting label transfers.
A query node with low match fraction is poorly matched to the reference
and its `transferred_label_accepted` will be `FALSE`.

``` r

diag <- somalign_diagnostics(fit)
diag$ot$match_fraction   # per query node
#> [1] 0.61832217 0.45163748 0.85627399 0.05511686 0.73599227 0.78775970
diag$ot$transport_mass   # total transported mass (< 1 = some mass discarded)
#> [1] 0.6046437
```

Use
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
to check that findings are stable across OT hyperparameters:

``` r

somalign_sensitivity_grid(
  query, reference,
  epsilon   = c(0.05, 0.1),
  rho_query = c(0.5, 1),
  rho_ref   = 1,
  solver    = "internal"
)
#>   epsilon rho_query rho_ref   solver transport_mass mean_match_fraction
#> 1    0.05       0.5       1 internal      0.5509337           0.5195532
#> 2    0.10       0.5       1 internal      0.5963177           0.5862370
#> 3    0.05       1.0       1 internal      0.6046437           0.5841837
#> 4    0.10       1.0       1 internal      0.6418082           0.6353367
#>   max_row_mass_error max_col_mass_error accepted_label_fraction
#> 1          0.2094605          0.1317681               0.6666667
#> 2          0.2089333          0.1228897               0.6666667
#> 3          0.1989228          0.1239883               1.0000000
#> 4          0.1973681          0.1173270               0.8333333
#>   outside_direct_fraction outside_corrected_fraction
#> 1               0.2631579                  0.2894737
#> 2               0.2631579                  0.2894737
#> 3               0.2631579                  0.1052632
#> 4               0.2631579                  0.1052632
```

## Practical checklist

- Use identical feature names and feature order.
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  reorders columns to `reference$features` and fails if any are missing.
- Train the query SOM in reference-scaled space as shown above.
- Treat `corrected_som_unit` and `transferred_label` as annotation
  support. Keep `old_som_unit`, `outside_reference_distance`, and
  `final_status` as the conservative classification.
- If you only have a saved query SOM codebook, you can pass a matrix as
  `som_query`, provided it is already in the reference-scaled feature
  space with columns named with `reference$features`.
