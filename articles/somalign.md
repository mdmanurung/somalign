# Quick start

`somalign` aligns a query self-organising map to a fixed reference SOM
using codebook-level unbalanced entropic optimal transport. The example
below trains a reference SOM on labelled old samples, then projects a
shifted query dataset into that reference.

## Quick start

``` r

library(kohonen)
library(somalign)

set.seed(1)
old <- rbind(
  matrix(rnorm(80, mean = -1), ncol = 4),
  matrix(rnorm(80, mean = 1), ncol = 4)
)
colnames(old) <- paste0("f", seq_len(ncol(old)))
labels <- rep(c("low", "high"), each = 20)

reference <- somalign_train_reference(
  old,
  labels = labels,
  grid = kohonen::somgrid(2, 2, "hexagonal"),
  rlen = 5
)

query <- old + 0.2
query_obj <- somalign_query(
  query,
  reference,
  grid = kohonen::somgrid(2, 2, "hexagonal"),
  rlen = 5
)

fit <- somalign_fit(query_obj, reference)
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.94); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
results <- somalign_results(fit)
results_with_meta <- somalign_results(
  fit,
  data = data.frame(batch = rep("query_1", nrow(query)))
)
```

## Current API surface

The user-facing workflow is built around a small set of exported
functions:

- Build a reference with
  [`somalign_train_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_train_reference.md),
  or wrap saved artifacts with
  [`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
  /
  [`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md).
- Prepare query data with
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).
  Query samples are always scaled with the reference center and scale.
- Align query and reference SOM nodes with
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  The default internal solver is used for both `solver = "internal"` and
  the compatibility alias `solver = "auto"`.
- Extract per-sample output with
  [`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md).
  Pass `somalign_results(fit, data = sample_metadata)` to append one
  metadata row per query sample.
- Inspect fit quality with
  [`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
  and tune OT parameters with
  [`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md).

## What to look at first

[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
returns one row per query sample. Start with these columns:

``` r

head(results[, c(
  "sample_id",
  "query_som_unit",
  "old_som_unit",
  "old_som_distance",
  "old_som_distance_threshold",
  "old_som_label",
  "outside_reference_distance",
  "final_status",
  "corrected_som_unit",
  "corrected_som_distance",
  "corrected_som_distance_threshold",
  "corrected_outside_reference_distance",
  "correction_norm",
  "transferred_label",
  "transferred_label_confidence",
  "transferred_label_accepted"
)])
#>   sample_id query_som_unit old_som_unit old_som_distance
#> 1         1              4            1        1.6019468
#> 2         2              1            4        0.8377855
#> 3         3              1            4        1.3090895
#> 4         4              3            4        1.9063036
#> 5         5              1            4        0.8813879
#> 6         6              1            4        0.9746628
#>   old_som_distance_threshold old_som_label outside_reference_distance
#> 1                   1.775302          high                      FALSE
#> 2                   1.919417           low                      FALSE
#> 3                   1.919417           low                      FALSE
#> 4                   1.919417           low                      FALSE
#> 5                   1.919417           low                      FALSE
#> 6                   1.919417           low                      FALSE
#>       final_status corrected_som_unit corrected_som_distance
#> 1 inside_reference                  3              1.1371940
#> 2 inside_reference                  4              1.1086155
#> 3 inside_reference                  4              1.5346506
#> 4 inside_reference                  4              1.5547127
#> 5 inside_reference                  4              0.9885232
#> 6 inside_reference                  4              0.5988922
#>   corrected_som_distance_threshold corrected_outside_reference_distance
#> 1                         1.532219                                FALSE
#> 2                         1.919417                                FALSE
#> 3                         1.919417                                FALSE
#> 4                         1.919417                                FALSE
#> 5                         1.919417                                FALSE
#> 6                         1.919417                                FALSE
#>   correction_norm transferred_label transferred_label_confidence
#> 1       1.3752534              <NA>                           NA
#> 2       0.6492682               low                    0.8433909
#> 3       0.6492682               low                    0.8433909
#> 4       0.8008908              <NA>                           NA
#> 5       0.6492682               low                    0.8433909
#> 6       0.6492682               low                    0.8433909
#>   transferred_label_accepted
#> 1                      FALSE
#> 2                       TRUE
#> 3                       TRUE
#> 4                      FALSE
#> 5                       TRUE
#> 6                       TRUE
```

## Interpreting the result

**Direct projection** (`old_som_unit`, `old_som_distance`,
`old_som_distance_threshold`, `outside_reference_distance`,
`final_status`, `old_som_label`, `old_som_label_confidence`) is the
primary result. Each query sample is assigned to its nearest reference
node by Euclidean distance in the shared feature space. No transport is
involved; the assignment is deterministic given the reference codebook.

**Corrected projection** (`corrected_som_unit`,
`corrected_som_distance`, `corrected_som_distance_threshold`,
`corrected_outside_reference_distance`, `correction_norm`) is auxiliary.
The OT plan shifts each query SOM node toward its matched reference
nodes, and query samples are then re-projected from the shifted
position. A large `correction_norm` points to a systematic displacement
between the query and reference populations — useful for visualisation
and triage, but not a replacement for the direct assignment.

**Transferred labels** (`transferred_label`,
`transferred_label_confidence`, `transferred_label_accepted`) are also
auxiliary, derived from the same OT correspondence. A label is accepted
only when the query node has sufficient transported mass and the
dominant matched reference label is confident enough. Always check
`transferred_label_accepted` before using a transferred label
downstream.

## Next steps

[`vignette("pretrained-old-and-new-soms", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/pretrained-old-and-new-soms.md)
covers the workflow for existing reference or query SOMs:
reference-scaled codebooks, node-level artifacts, diagnostics, and
hyperparameter tuning.

[`vignette("anchor-samples", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/anchor-samples.md)
covers
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md):
supplying remeasured QC samples as anchor pairs, tuning `rho_anchor`,
and inspecting anchor node coverage.

[`vignette("algorithm", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/algorithm.md)
walks through each pipeline stage — direct projection, OT
correspondence, correction vectors, label transfer — and explains how
they produce the output columns.

For fitted objects,
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
reports solver convergence, OT mass behaviour, node-level match
fractions, and direct/corrected outside reference fractions.
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
sweeps OT hyperparameters to check whether corrected projection and
label transfer are stable. For large query matrices,
`somalign_fit(chunk_size = ...)` controls how many samples are projected
at once during nearest-reference-node searches.

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
    #>  [1] cli_3.6.6           knitr_1.51          rlang_1.3.0        
    #>  [4] xfun_0.60           otel_0.2.0          textshaping_1.0.5  
    #>  [7] jsonlite_2.0.0      htmltools_0.5.9     ragg_1.5.2         
    #> [10] sass_0.4.10         rmarkdown_2.31      evaluate_1.0.5     
    #> [13] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
    #> [16] lifecycle_1.0.5     bookdown_0.47       BiocManager_1.30.27
    #> [19] compiler_4.6.1      fs_2.1.0            htmlwidgets_1.6.4  
    #> [22] Rcpp_1.1.2          systemfonts_1.3.2   digest_0.6.39      
    #> [25] R6_2.6.1            bslib_0.11.0        tools_4.6.1        
    #> [28] pkgdown_2.2.1       cachem_1.1.0        desc_1.4.3
