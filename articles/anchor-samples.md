# Anchor-regularized alignment with remeasured samples

In many longitudinal and multi-site studies, a small set of
quality-control samples – control cell lines, reference bead
populations, or proficiency panel specimens – gets run alongside every
batch. When those same samples appear in both the old batch (from which
the reference SOM was trained) and the new batch (from which the query
SOM was trained), they carry direct information about the per-node batch
displacement: wherever an anchor sample lands in the reference codebook,
the corresponding query measurement tells you how far that region has
shifted.

[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
folds that information into the optimal transport solve. For each anchor
pair, the old measurement is projected onto the query SOM and the new
measurement onto the reference SOM. (The query SOM was trained on
new-batch data, so projecting the old-batch anchor onto it identifies
which query node that anchor occupied before the batch shift; projecting
the new-batch anchor onto the reference SOM identifies the corresponding
reference node after the shift.) This produces a node-pair count matrix
– how many anchor pairs link query node *k* to reference node *l*. Pairs
with anchor support get reduced transport cost, so the OT plan
preferentially routes mass through those node combinations. Uncovered
nodes remain unaffected; their cost is identical to the standard
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
solve.

## Setup

``` r

library(kohonen)
library(somalign)

set.seed(42)
p <- 5L

# Old batch: two well-separated populations
old_data <- rbind(
  matrix(rnorm(60 * p, mean = -2, sd = 0.6), ncol = p),
  matrix(rnorm(60 * p, mean =  2, sd = 0.6), ncol = p)
)
colnames(old_data) <- paste0("CD", seq_len(p))

# Train reference SOM and build reference object in one call
reference <- somalign_train_reference(
  old_data,
  grid = kohonen::somgrid(4, 4, "hexagonal"),
  rlen = 50
)
reference
#> <somalign_reference>
#>   features: 5
#>   reference nodes: 16
#>   labelled nodes: 0
```

The new batch carries a modest uniform batch offset – realistic for
reagent lot changes or inter-instrument calibration drift.

``` r

shift    <- rep(0.2, p)
new_data <- old_data + matrix(shift, nrow = nrow(old_data), ncol = p, byrow = TRUE)

query <- somalign_query(
  new_data, reference,
  grid = kohonen::somgrid(4, 4, "hexagonal"),
  rlen = 50
)
query
#> <somalign_query>
#>   samples: 120
#>   features: 5
#>   query nodes: 16
```

## Identifying anchor samples

Suppose 25 samples from a quality-control pool – the same biological
material measured under both batch conditions – are spread across both
phenotypic populations.

``` r

# 12 anchors from the first population, 13 from the second
anc_idx    <- c(seq_len(12L), 60L + seq_len(13L))
anchor_old <- old_data[anc_idx, ]   # QC pool in the old batch
anchor_new <- new_data[anc_idx, ]   # Same QC pool in the new batch

cat("Anchor pairs:", nrow(anchor_old), "\n")
#> Anchor pairs: 25
cat("Features:    ", ncol(anchor_old), "\n")
#> Features:     5
```

## Comparing standard and anchored fits

``` r

fit_plain <- somalign_fit(query, reference)
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.49); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
fit_anc   <- somalign_fit_anchored(query, reference,
                                    anchor_old = anchor_old,
                                    anchor_new = anchor_new,
                                    rho_anchor = 1.5)
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.52); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
```

Printing either object gives a concise summary:

``` r

fit_plain
#> <somalign_fit>
#>   solver: internal
#>   query nodes: 16
#>   reference nodes: 16
#>   transport mass: 2.094
fit_anc
#> <somalign_anchored_fit>
#>   solver: internal
#>   query nodes: 16
#>   reference nodes: 16
#>   transport mass: 2.114
#>   anchors: 25 (75% node coverage)
```

The anchored fit reports how many query nodes had at least one anchor
pair and the fraction of the codebook they cover:

``` r

fit_anc$anchors
#> $n_anchors
#> [1] 25
#> 
#> $rho_anchor
#> [1] 1.5
#> 
#> $nodes_covered
#> [1] 12
#> 
#> $coverage_fraction
#> [1] 0.75
```

`nodes_covered` counts distinct query SOM nodes onto which at least one
anchor sample’s old-batch measurement was projected. `coverage_fraction`
is that count divided by the total number of query nodes. Nodes with no
anchor coverage are solved with unchanged costs, so they behave exactly
as they would under
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## Effect on node correction vectors

Anchor regularization shifts mass toward node pairs that anchors
confirm. A practical way to see this is to compare the per-node
correction magnitudes between the two fits. Each correction vector
points from a query SOM node toward its barycentric target in the
reference, so nodes whose OT correspondence changed will show different
norms:

``` r

plain_norms <- sqrt(rowSums(fit_plain$node_shifts^2))
anc_norms   <- sqrt(rowSums(fit_anc$node_shifts^2))

cat("Mean correction norm -- plain:    ", round(mean(plain_norms), 4), "\n")
#> Mean correction norm -- plain:     0.3801
cat("Mean correction norm -- anchored: ", round(mean(anc_norms),   4), "\n")
#> Mean correction norm -- anchored:  0.3734
cat("Max  correction norm -- plain:    ", round(max(plain_norms),  4), "\n")
#> Max  correction norm -- plain:     0.6273
cat("Max  correction norm -- anchored: ", round(max(anc_norms),    4), "\n")
#> Max  correction norm -- anchored:  0.623
```

When anchors consistently confirm the same direction of displacement,
anchor-regularized fits produce more coherent correction vectors: the OT
plan concentrates mass on anchor-supported routes, reducing the
influence of diffuse entropic transport on uncovered node pairs.

## Downstream results

Both fit objects expose the same downstream interface.
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
returns the standard projection data frame for all query samples:

``` r

res_plain <- somalign_results(fit_plain)
res_anc   <- somalign_results(fit_anc)

# Fraction of samples landing outside reference distance thresholds
cat("Outside fraction (plain):    ",
    mean(res_plain$final_status == "outside_reference"), "\n")
#> Outside fraction (plain):     0.2333333
cat("Outside fraction (anchored): ",
    mean(res_anc$final_status   == "outside_reference"), "\n")
#> Outside fraction (anchored):  0.2333333
```

[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
works identically on anchored fits:

``` r

diag <- somalign_diagnostics(fit_anc)
cat("Solver converged:", diag$solver$converged, "\n")
#> Solver converged: TRUE
cat("Transport mass:  ", round(diag$ot$transport_mass, 4), "\n")
#> Transport mass:   2.1137
```

## Tuning `rho_anchor`

`rho_anchor` controls how much the anchor bonus reduces the normalized
OT cost. At `rho_anchor = 0`, the bonus vanishes and the result is
identical to
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
As `rho_anchor` grows, anchor-supported pairs become progressively
cheaper. Once the bonus exceeds the normalized cost for a pair, that
pair’s effective cost is clamped to zero – the plan then spreads mass
among all zero-cost pairs according to entropic regularization rather
than differential anchor support, so very large values do not continue
to concentrate mass.

A practical upper bound is `rho_anchor * max_A / n_anchors <= 1`, where
`max_A` is the highest anchor count for any single node pair. Values in
the range 0.5–2 are a reasonable starting point. When coverage is sparse
(\< 30% of nodes) or anchor pairs are unevenly distributed, a smaller
value avoids over-weighting a few node pairs.

``` r

rhos <- c(0.5, 1, 1.5, 2, 4)
sweep_results <- lapply(rhos, function(rho) {
  f <- somalign_fit_anchored(query, reference,
                              anchor_old = anchor_old,
                              anchor_new = anchor_new,
                              rho_anchor = rho)
  data.frame(
    rho_anchor     = rho,
    transport_mass = round(f$diagnostics$ot$transport_mass, 5),
    mean_corr_norm = round(mean(sqrt(rowSums(f$node_shifts^2))), 5)
  )
})
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.49); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.50); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.52); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.53); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 14 query node(s) have match_mass_ratio > 1 (max 4.60); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
do.call(rbind, sweep_results)
#>   rho_anchor transport_mass mean_corr_norm
#> 1        0.5        2.10163        0.37732
#> 2        1.0        2.10824        0.37503
#> 3        1.5        2.11371        0.37341
#> 4        2.0        2.11936        0.37181
#> 5        4.0        2.13675        0.36734
```

For a more targeted diagnostic, compare the
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
node-level `match_fraction` and `correction_norm` across values, and use
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
to check whether corrected node assignments shift in a direction
consistent with your biological expectations.

## When to use `somalign_fit_anchored()`

Use this function when:

- You have samples that were physically run in both the old and new
  batches (QC pools, reference standards, proficiency panels).
- You want the OT plan to respect known per-node correspondence, rather
  than relying entirely on node-mass distributions.

Anchor pairs need not cover the entire dataset – even sparse coverage
meaningfully constrains the OT plan for covered node pairs, while
uncovered nodes fall back to the standard OT objective.

If no remeasured samples are available,
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
is the correct choice. There is no benefit to calling the anchored
variant with fabricated or imputed anchor pairs.

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
