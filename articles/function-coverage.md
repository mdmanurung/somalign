# API smoke test — all public functions

This vignette exercises every public function in the package so that
build-time vignette rendering catches breaking changes or silent
regressions. It is not a tutorial — see the other vignettes for
annotated workflows. Each section calls the function, checks the return
type, and prints a key field to confirm the object was constructed
correctly.

``` r

library(kohonen)
library(somalign)

set.seed(42)
p <- 4L

old_data <- rbind(
  matrix(rnorm(30 * p, mean = -2, sd = 0.5), ncol = p),
  matrix(rnorm(30 * p, mean =  2, sd = 0.5), ncol = p)
)
colnames(old_data) <- paste0("M", seq_len(p))
labels_old <- rep(c("pop_A", "pop_B"), each = 30)

new_data <- old_data + 0.3
grid_small <- kohonen::somgrid(3, 3, "hexagonal")
```

## `somalign_train_reference`

Trains a SOM on old-batch data and returns a `somalign_reference`.

``` r

ref <- somalign_train_reference(
  old_data, labels = labels_old,
  grid = grid_small, rlen = 20
)
stopifnot(inherits(ref, "somalign_reference"))
ref
#> <somalign_reference>
#>   features: 4
#>   reference nodes: 9
#>   labelled nodes: 7
```

## `somalign_reference`

Builds a reference from an *existing* `kohonen` SOM object rather than
training a new one. This path is used when the reference SOM was trained
externally or loaded from a saved file.

``` r

scaled_old <- scale(old_data)
som_obj    <- kohonen::som(scaled_old, grid = grid_small, rlen = 20)
ref2 <- somalign_reference(
  som_ref = som_obj,
  data    = old_data,
  labels  = labels_old,
  codebook_space = "reference_scaled"
)
stopifnot(inherits(ref2, "somalign_reference"))
cat("n_samples:", ref2$n_samples, "\n")
#> n_samples: 60
```

## `somalign_reference_from_nodes`

Reconstructs a reference from stored node-level artifacts. This is
useful when loading a pre-computed reference from disk without
re-running the SOM.

``` r

ref3 <- somalign_reference_from_nodes(
  codebook                  = ref$codebook,
  features                  = ref$features,
  center                    = ref$center,
  scale                     = ref$scale,
  node_masses               = ref$node_masses,
  label_prob                = ref$label_prob,
  distance_quantiles        = ref$distance_quantiles,
  global_distance_quantiles = ref$global_distance_quantiles
)
stopifnot(inherits(ref3, "somalign_reference"))
stopifnot(identical(ref3$features, ref$features))
cat("Nodes:", nrow(ref3$codebook), "\n")
#> Nodes: 9
```

## `somalign_reference_from_som`

Builds a reference from a
[`kohonen::supersom()`](https://rdrr.io/pkg/kohonen/man/supersom.html)
or [`kohonen::xyf()`](https://rdrr.io/pkg/kohonen/man/supersom.html) SOM
that already has a label code layer. Requires `keep.data = TRUE` so that
cell-level unit assignments and raw data are accessible for per-cell
distance computation. With `labels = "codebook"`, the second code book
populates `label_prob`, enabling `label_guided` alignment.

``` r

Y_old   <- model.matrix(~ labels_old - 1)
colnames(Y_old) <- levels(factor(labels_old))
old_scaled <- scale(old_data, center = ref$center, scale = ref$scale)

ss_ref <- kohonen::supersom(
  list(old_scaled, Y_old),
  grid      = grid_small,
  rlen      = 20,
  keep.data = TRUE
)
ref_ss <- somalign_reference_from_som(
  ss_ref,
  center         = ref$center,
  scale          = ref$scale,
  codebook_space = "reference_scaled",
  labels         = "codebook"
)
stopifnot(inherits(ref_ss, "somalign_reference"))
stopifnot(!is.null(ref_ss$label_prob))
cat("Nodes:", nrow(ref_ss$codebook), "  label_prob rows:", nrow(ref_ss$label_prob), "\n")
#> Nodes: 9   label_prob rows: 9
```

## `somalign_query`

Trains a query SOM on new-batch data and projects it against the
reference.

``` r

qry <- somalign_query(new_data, ref, grid = grid_small, rlen = 20)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
stopifnot(inherits(qry, "somalign_query"))
qry
#> <somalign_query>
#>   samples: 60
#>   features: 4
#>   query nodes: 9
```

## `somalign_query_from_som`

Builds a query object from a pre-trained `kohonen` SOM without
re-computing nearest-code assignments. Reuses `som$unit.classif`
directly, so `sample_distance` is set to NA; all downstream alignment
functions still work.

``` r

new_scaled <- scale(new_data, center = ref$center, scale = ref$scale)
som_qry    <- kohonen::som(new_scaled, grid = grid_small, rlen = 20)

qry_reuse <- somalign_query_from_som(
  som_qry, data = new_data, reference = ref
)
stopifnot(inherits(qry_reuse, "somalign_query"))
cat("Sample units range:", range(qry_reuse$sample_unit), "\n")
#> Sample units range: 1 9
cat("sample_distance NA:", all(is.na(qry_reuse$sample_distance)), "\n")
#> sample_distance NA: TRUE
```

## `somalign_normalize`

Pre-aligns query data to the reference mean (and optionally scale) in
z-scored space. Returns a matrix in the original (raw) units, suitable
as `data` for
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

``` r

new_normed <- somalign_normalize(new_data, reference = ref, method = "mean")
stopifnot(is.matrix(new_normed))
stopifnot(identical(dim(new_normed), dim(new_data)))
cat("Max absolute shift removed:", round(max(abs(new_data - new_normed)), 4), "\n")
#> Max absolute shift removed: 0.3
```

## `somalign_quantile_normalize`

Divides each feature by its upper quantile in raw space, mapping the
bulk signal to approximately \[0, 1\]. For features with very different
dynamic ranges across batches.

``` r

new_qnormed <- somalign_quantile_normalize(new_data, reference = ref, probs = 0.999)
stopifnot(is.matrix(new_qnormed))
stopifnot(identical(dim(new_qnormed), dim(new_data)))
cat("Max value before:", round(max(new_data), 4), "\n")
#> Max value before: 3.3298
cat("Max value after: ", round(max(new_qnormed), 4), "\n")
#> Max value after:  1.005
```

## `somalign_fit`

Standard (unanchored) unbalanced entropic OT alignment.

``` r

fit <- somalign_fit(qry, ref)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(fit, "somalign_fit"))
fit
#> <somalign_fit>
#>   solver: internal
#>   query nodes: 9
#>   reference nodes: 9
#>   transport mass: 1.104
```

## `somalign_fit_anchored`

Anchor-regularized alignment using remeasured QC samples.

``` r

anc_idx    <- c(1:10, 31:40)
anchor_old <- old_data[anc_idx, ]
anchor_new <- new_data[anc_idx, ]

fit_anc <- somalign_fit_anchored(
  qry, ref,
  anchor_old = anchor_old,
  anchor_new = anchor_new,
  rho_anchor = 1.0
)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.14); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(fit_anc, "somalign_anchored_fit"))
stopifnot(inherits(fit_anc, "somalign_fit"))
cat("n_anchors:        ", fit_anc$anchors$n_anchors, "\n")
#> n_anchors:         20
cat("coverage_fraction:", round(fit_anc$anchors$coverage_fraction, 3), "\n")
#> coverage_fraction: 0.667
```

`rho_anchor = 0` short-circuits to the plain-OT path and emits a
message:

``` r

fit_anc0 <- somalign_fit_anchored(
  qry, ref,
  anchor_old = anchor_old,
  anchor_new = anchor_new,
  rho_anchor = 0
)
#> `rho_anchor = 0`: the anchor cost bonus is inactive. Use `somalign_fit()` for equivalent results.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(fit_anc0, "somalign_anchored_fit"))
stopifnot(fit_anc0$anchors$rho_anchor == 0)
```

`correction = "subspace"` estimates batch directions from the anchor
displacement SVD and restricts node shifts to that subspace, preserving
orthogonal biology:

``` r

fit_sub <- somalign_fit_anchored(
  qry, ref,
  anchor_old = anchor_old,
  anchor_new = anchor_new,
  correction = "subspace"
)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(fit_sub, "somalign_anchored_fit"))
stopifnot(!is.null(fit_sub$anchors$batch_subspace))
cat("Correction mode:    ", fit_sub$anchors$correction, "\n")
#> Correction mode:     subspace
cat("Batch subspace rank:", fit_sub$anchors$batch_subspace$rank, "\n")
#> Batch subspace rank: 1
cat("Variance explained: ",
    round(fit_sub$anchors$batch_subspace$variance_explained, 3), "\n")
#> Variance explained:  1
```

## `somalign_fit_two_pass`

Two-pass alignment: a first OT pass at `epsilon_global` estimates a
mass-weighted global shift, subtracts it from the query codebook, then a
second pass at `epsilon_local` captures residual per-node offsets.
Exposes `$two_pass` in addition to the standard `somalign_fit`
interface.

``` r

fit2 <- somalign_fit_two_pass(qry, ref, epsilon_global = 0.3, epsilon_local = 0.1)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.49); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.14); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(fit2, "somalign_fit"))
stopifnot(!is.null(fit2$two_pass))
stopifnot(all(c("global_shift", "global_shift_norm",
                "epsilon_global", "epsilon_local",
                "batch_subspace") %in% names(fit2$two_pass)))
cat("Global shift norm:  ", round(fit2$two_pass$global_shift_norm, 4), "\n")
#> Global shift norm:   0.2949
cat("Batch subspace rank:", fit2$two_pass$batch_subspace$rank, "\n")
#> Batch subspace rank: 2
```

[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
and
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
accept the two-pass fit unchanged:

``` r

res2  <- somalign_results(fit2)
diag2 <- somalign_diagnostics(fit2)
stopifnot(is.data.frame(res2))
stopifnot(all(c("solver", "ot", "nodes", "projection") %in% names(diag2)))
cat("Converged:", diag2$solver$converged, "\n")
#> Converged: TRUE
```

## `somalign_results`

Extracts the per-sample projection data frame from a fit. Works on both
plain and anchored fits.

``` r

res     <- somalign_results(fit)
res_anc <- somalign_results(fit_anc)

stopifnot(is.data.frame(res))
stopifnot(is.data.frame(res_anc))
stopifnot("corrected_som_unit" %in% names(res))
stopifnot("final_status" %in% names(res))
cat("Rows:", nrow(res), "  Columns:", ncol(res), "\n")
#> Rows: 60   Columns: 17

res_aug <- somalign_results(fit, data = new_data)
stopifnot(ncol(res_aug) > ncol(res))
```

## `somalign_diagnostics`

Returns a named list of solver, OT, node, and projection diagnostics.

``` r

diag <- somalign_diagnostics(fit)
stopifnot(is.list(diag))
stopifnot(all(c("solver", "ot", "nodes", "projection") %in% names(diag)))
cat("Converged:     ", diag$solver$converged, "\n")
#> Converged:      TRUE
cat("Transport mass:", round(diag$ot$transport_mass, 4), "\n")
#> Transport mass: 1.1036

# Same interface on anchored fit
diag_anc <- somalign_diagnostics(fit_anc)
stopifnot(!is.null(diag_anc$solver$converged))
```

## `somalign_sensitivity_grid`

Sweeps `epsilon` × `rho_query` × `rho_ref` and returns a summary data
frame.

``` r

sg <- somalign_sensitivity_grid(
  qry, ref,
  epsilon   = c(0.1, 0.3, 0.7),
  rho_query = c(0.5, 1.0),
  rho_ref   = 1.0
)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.17); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.68); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 2.83); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.13); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.49); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 2.32); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(is.data.frame(sg))
stopifnot(nrow(sg) == 6L)
stopifnot("transport_mass" %in% names(sg))
sg[, c("epsilon", "rho_query", "transport_mass")]
#>   epsilon rho_query transport_mass
#> 1     0.1       0.5       1.138275
#> 2     0.3       0.5       1.501523
#> 3     0.7       0.5       2.268484
#> 4     0.1       1.0       1.103620
#> 5     0.3       1.0       1.373678
#> 6     0.7       1.0       1.945012
```

## `somalign_som_stability`

Measures how much OT statistics vary with query SOM random seed by
re-training the query SOM from scratch for each seed.

``` r

stab <- somalign_som_stability(
  new_data, ref,
  som_seeds = 1:3,
  grid = grid_small,
  rlen = 20
)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.27); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.22); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.14); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(is.data.frame(stab))
stopifnot(nrow(stab) == 3L)
stopifnot("som_seed" %in% names(stab))
stab[, c("som_seed", "transport_mass", "converged")]
#>   som_seed transport_mass converged
#> 1        1       1.101402      TRUE
#> 2        2       1.111093      TRUE
#> 3        3       1.110050      TRUE
```

## Print methods

All four main classes have print methods that should return the object
invisibly.

``` r

invisible(print(ref))
#> <somalign_reference>
#>   features: 4
#>   reference nodes: 9
#>   labelled nodes: 7
invisible(print(qry))
#> <somalign_query>
#>   samples: 60
#>   features: 4
#>   query nodes: 9
invisible(print(fit))
#> <somalign_fit>
#>   solver: internal
#>   query nodes: 9
#>   reference nodes: 9
#>   transport mass: 1.104
invisible(print(fit_anc))
#> <somalign_anchored_fit>
#>   solver: internal
#>   query nodes: 9
#>   reference nodes: 9
#>   transport mass: 1.123
#>   anchors: 20 (66.7% node coverage)
```

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
