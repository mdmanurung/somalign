# Advanced workflow: existing SOMs

Many analysis pipelines train a reference SOM once and save it, then
process new batches against that fixed reference at a later stage.
`somalign` supports this: you can supply an existing reference SOM, a
separately trained query SOM, or pre-computed node-level artifacts
rather than letting `somalign` train everything from scratch.

## Coordinate-space requirements

The single constraint that ties the whole approach together is that both
the reference and query codebooks must live in the same feature
coordinate system: z-scored with the **old/reference** data’s mean and
standard deviation. When the reference SOM was trained on old data
scaled to $`\mu_\text{ref}`$ and $`\sigma_\text{ref}`$, the query SOM
must be trained on new samples transformed by those same parameters —
not by new-data z-scores, and not in raw units.

``` r

new_scaled_for_som <- sweep(
  sweep(new_data[, reference$features], 2, reference$center, "-"),
  2,
  reference$scale,
  "/"
)
new_som <- kohonen::som(new_scaled_for_som, grid = ...)
```

[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
still receives the original (unscaled) new data matrix; it applies the
reference scaling internally for the per-sample projection step. When
you pass a saved query SOM, set `codebook_space` to match the saved
codebook. The default is `"reference_scaled"`, which is correct for the
query SOM trained above. Use `"raw"` only when the saved query SOM
codebook is in raw feature units:

``` r

query <- somalign_query(
  new_data,
  reference,
  som_query = new_som_trained_on_raw_data,
  codebook_space = "raw"
)
```

## Building references from saved artifacts

[`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
is the preferred route when you have the old SOM object and the old
sample matrix. It recomputes node masses, per-node label probabilities,
and reference-distance thresholds from the data:

``` r

reference <- somalign_reference(
  old_som,
  old_data,
  codebook_space = "reference_scaled"
)
```

Set `codebook_space = "raw"` only when the saved SOM was trained on raw
(unscaled) feature values;
[`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
will apply the reference scaling to the codebook before use.

[`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md)
covers the case where neither the SOM object nor the original samples
are available, but node-level summaries were archived:

``` r

reference <- somalign_reference_from_nodes(
  codebook = old_codebook,
  features = feature_names,
  center = old_center,
  scale = old_scale,
  node_masses = old_node_masses,
  label_prob = old_label_prob,
  distance_quantiles = old_distance_quantiles
)
```

Without saved distance quantiles, the outside-reference distance flag
will be disabled; without saved label probabilities, label transfer will
be disabled. Direct node assignment continues to work in both cases.

## Full existing-SOM example

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
  som_query = new_som,
  codebook_space = "reference_scaled"
)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.

fit <- somalign_fit(query, reference)
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
sample_metadata <- data.frame(batch = rep("new_batch", nrow(new_data)))
results <- somalign_results(fit, data = sample_metadata)

head(results[, c(
  "sample_id",
  "query_som_unit",
  "old_som_unit",
  "old_som_label",
  "old_som_distance_threshold",
  "outside_reference_distance",
  "final_status",
  "corrected_som_unit",
  "corrected_som_distance_threshold",
  "corrected_outside_reference_distance",
  "correction_norm",
  "transferred_label",
  "transferred_label_confidence",
  "transferred_label_accepted"
)])
#>   sample_id query_som_unit old_som_unit old_som_label
#> 1         1              3            1       old_low
#> 2         2              6            5       old_low
#> 3         3              2            1       old_low
#> 4         4              3            1       old_low
#> 5         5              3            4       old_low
#> 6         6              3            1       old_low
#>   old_som_distance_threshold outside_reference_distance     final_status
#> 1                   2.483497                      FALSE inside_reference
#> 2                   1.255846                      FALSE inside_reference
#> 3                   2.483497                      FALSE inside_reference
#> 4                   2.483497                      FALSE inside_reference
#> 5                   2.023005                      FALSE inside_reference
#> 6                   2.483497                      FALSE inside_reference
#>   corrected_som_unit corrected_som_distance_threshold
#> 1                  1                         2.483497
#> 2                  5                         1.255846
#> 3                  5                         1.255846
#> 4                  1                         2.483497
#> 5                  4                         2.023005
#> 6                  1                         2.483497
#>   corrected_outside_reference_distance correction_norm transferred_label
#> 1                                FALSE       0.3303817           old_low
#> 2                                 TRUE       0.5424595              <NA>
#> 3                                 TRUE       1.1510211              <NA>
#> 4                                FALSE       0.3303817           old_low
#> 5                                FALSE       0.3303817           old_low
#> 6                                FALSE       0.3303817           old_low
#>   transferred_label_confidence transferred_label_accepted
#> 1                    0.9987403                       TRUE
#> 2                           NA                      FALSE
#> 3                           NA                      FALSE
#> 4                    0.9987403                       TRUE
#> 5                    0.9987403                       TRUE
#> 6                    0.9987403                       TRUE
```

## Quality control and tuning

Before acting on label transfer or corrected node assignments, inspect
the solver and OT diagnostics:

``` r

diagnostics <- somalign_diagnostics(fit)

diagnostics$solver[c(
  "used", "converged", "iterations", "final_delta", "cost_scale"
)]
#> $used
#> [1] "internal"
#> 
#> $converged
#> [1] TRUE
#> 
#> $iterations
#> [1] 80
#> 
#> $final_delta
#> [1] 9.047741e-08
#> 
#> $cost_scale
#> [1] 2.161104
diagnostics$ot$match_fraction       # clipped to at most 1
#> [1] 0.8271323 0.8592818 1.0000000 0.2727043 0.9308854 1.0000000
diagnostics$ot$match_mass_ratio     # raw transported/query mass ratio
#> [1] 0.8271323 0.8592818 1.0886089 0.2727043 0.9308854 1.1254354
diagnostics$ot$transport_mass       # total transported mass
#> [1] 0.8325581
diagnostics$ot$max_row_mass_error   # query-side marginal deviation
#> [1] 0.1531149
diagnostics$projection              # direct/corrected outside fractions
#> $outside_direct_fraction
#> [1] 0.2631579
#> 
#> $outside_corrected_fraction
#> [1] 0.1315789
```

`diagnostics$ot$match_fraction` is the most actionable single number: a
query node with low match fraction transported little mass to any
reference node, which typically means a novel population, a
codebook-space mismatch, or both. Labels transferred from such nodes
should be treated as provisional. `diagnostics$ot$match_mass_ratio`
stores the unclipped transported-mass ratio; values above 1 can occur in
unbalanced OT and are not errors.

`diagnostics$ot$transport_mass` tells you how much of the total query
mass the OT plan moved to reference nodes (as opposed to absorbing it
via the KL margin penalties). `diagnostics$ot$max_row_mass_error`
captures the worst per-node marginal deviation; a large value means the
unbalanced fit allowed substantial mass destruction, which is
appropriate for truly novel query content but should be visible in QC
records. The `diagnostics$projection` fractions report how many samples
are outside the saved reference distance thresholds before and after the
OT-derived correction.

Use
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
to confirm that auxiliary outputs hold across OT hyperparameters before
drawing conclusions:

``` r

sensitivity <- somalign_sensitivity_grid(
  query, reference,
  epsilon   = c(0.05, 0.1),
  rho_query = c(0.5, 1),
  rho_ref   = 1
)
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.07); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 1 query node(s) have match_mass_ratio > 1 (max 1.05); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
sensitivity[, c(
  "epsilon",
  "rho_query",
  "rho_ref",
  "solver",
  "transport_mass",
  "mean_match_fraction",
  "max_row_mass_error",
  "accepted_label_fraction",
  "outside_direct_fraction",
  "outside_corrected_fraction"
)]
#>   epsilon rho_query rho_ref   solver transport_mass mean_match_fraction
#> 1    0.05       0.5       1 internal      0.7567440           0.7269262
#> 2    0.10       0.5       1 internal      0.8146693           0.7742114
#> 3    0.05       1.0       1 internal      0.7866239           0.7705905
#> 4    0.10       1.0       1 internal      0.8325581           0.8150006
#>   max_row_mass_error accepted_label_fraction outside_direct_fraction
#> 1          0.1925167               0.8333333               0.2631579
#> 2          0.1885548               0.6666667               0.2631579
#> 3          0.1576913               0.6666667               0.2631579
#> 4          0.1531149               0.6666667               0.2631579
#>   outside_corrected_fraction
#> 1                 0.10526316
#> 2                 0.10526316
#> 3                 0.07894737
#> 4                 0.13157895
```

Strong sensitivity to `epsilon` or `rho_query` is a signal that the
OT-derived outputs are not reliable for this dataset. In that case, rely
on direct projection and do not carry transferred labels or corrected
assignments into downstream analyses.

## Before relying on results

Feature names and order must match exactly —
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
reorders columns to `reference$features` and will error if any are
missing. Both codebooks must live in the same feature space; the query
SOM must be trained with `reference$center` and `reference$scale` unless
you supply a raw saved SOM with `codebook_space = "raw"`.

`old_som_unit`, `outside_reference_distance`, `final_status`, and
`old_som_label` are the primary result; corrected projection and
transferred labels are auxiliary. Use
`somalign_results(fit, data = sample_metadata)` to attach one metadata
row per query sample. Examine
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
before trusting label transfer, and run
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
whenever conclusions rest on corrected nodes or transferred labels.

A raw matrix of codebook vectors can be passed as `som_query` directly
if it is already in reference-scaled space with columns named to match
`reference$features`. For large query matrices,
`somalign_fit(chunk_size = ...)` controls how many samples are projected
at once.

## Session info

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] somalign_0.99.1  kohonen_3.0.13   BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
    #>  [4] bookdown_0.47       fastmap_1.2.0       xfun_0.60          
    #>  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
    #> [10] rmarkdown_2.31      lifecycle_1.0.5     cli_3.6.6          
    #> [13] sass_0.4.10         pkgdown_2.2.1       jquerylib_0.1.4    
    #> [16] compiler_4.6.1      tools_4.6.1         bslib_0.11.0       
    #> [19] evaluate_1.0.5      Rcpp_1.1.2          yaml_2.3.12        
    #> [22] BiocManager_1.30.27 otel_0.2.0          jsonlite_2.0.0     
    #> [25] rlang_1.3.0         fs_2.1.0            htmlwidgets_1.6.4
