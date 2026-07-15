# Two-pass global and local alignment

Batch effects in cytometry often contain a dominant global component —
reagent lot changes, instrument calibration drift, or environmental
factors that shift all populations in roughly the same direction by
roughly the same amount — layered on top of smaller population-specific
displacements. A single OT pass at a small epsilon captures fine
structure well but can be destabilised when the global offset dwarfs the
typical inter-node distance. A single pass at large epsilon smooths away
that instability at the cost of resolution.

[`somalign_fit_two_pass()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_two_pass.md)
resolves this explicitly: a first OT pass at a larger `epsilon_global`
estimates the mass-weighted mean node displacement (the *global shift*),
subtracts it from the query codebook, and then a second pass at a finer
`epsilon_local` captures residual population-specific offsets. The final
per-node correction is the sum of both passes, so the total displacement
for each sample is equivalent to a direct single-pass alignment at
`epsilon_local` — but the global subtraction makes the second-pass OT
problem better conditioned.

## Setup

``` r

library(kohonen)
library(somalign)

set.seed(42)
p <- 4L

old_data <- rbind(
  matrix(rnorm(60 * p, mean = -2, sd = 0.5), ncol = p),
  matrix(rnorm(60 * p, mean =  2, sd = 0.5), ncol = p)
)
colnames(old_data) <- paste0("CD", seq_len(p))
old_labels <- rep(c("pop_A", "pop_B"), each = 60L)

reference <- somalign_train_reference(
  old_data,
  labels = old_labels,
  grid   = kohonen::somgrid(3, 3, "hexagonal"),
  rlen   = 30
)
reference
#> <somalign_reference>
#>   features: 4
#>   reference nodes: 9
#>   labelled nodes: 7
```

The new batch carries a clear global offset plus smaller
within-population noise:

``` r

global_shift <- c(0.4, 0.3, 0.35, 0.25)

set.seed(43)
local_noise <- matrix(rnorm(nrow(old_data) * p, 0, 0.1), ncol = p)
new_data    <- old_data +
  matrix(global_shift, nrow(old_data), p, byrow = TRUE) + local_noise

query <- somalign_query(
  new_data, reference,
  grid = kohonen::somgrid(3, 3, "hexagonal"),
  rlen = 30
)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
query
#> <somalign_query>
#>   samples: 120
#>   features: 4
#>   query nodes: 9
```

## Two-pass alignment

``` r

fit2 <- somalign_fit_two_pass(
  query, reference,
  epsilon_global = 0.3,
  epsilon_local  = 0.1
)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.41); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.08); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
fit2
#> <somalign_fit>
#>   solver: internal
#>   query nodes: 9
#>   reference nodes: 9
#>   transport mass: 1.053
```

The `$two_pass` slot records the estimated global shift and both epsilon
values:

``` r

cat("Estimated global shift (correction direction, query -> reference):\n")
#> Estimated global shift (correction direction, query -> reference):
print(round(fit2$two_pass$global_shift, 3))
#>    CD1    CD2    CD3    CD4 
#> -0.165 -0.132 -0.153 -0.094
cat("\nTrue batch offset (positive = new batch above old):\n")
#> 
#> True batch offset (positive = new batch above old):
print(round(global_shift, 3))
#> [1] 0.40 0.30 0.35 0.25
cat("\nGlobal shift magnitude:", round(fit2$two_pass$global_shift_norm, 3), "\n")
#> 
#> Global shift magnitude: 0.277
```

`global_shift` points in the *correction* direction — from query
codebook toward the reference — so its sign is opposite to the true
batch offset. Because the node shifts are the mass-weighted mean of
first-pass OT transport vectors, the magnitude is in reference-scaled
units, not raw units, and does not directly equal the raw feature
offset. What matters in practice is the direction: each feature’s sign
should oppose the known batch offset, and features with larger
displacements should show larger absolute values. The norm
`global_shift_norm` gives an overall scalar measure of the global batch
displacement that can be compared across batches or experiments.

## Batch subspace diagnostic

The `$two_pass$batch_subspace` element summarises the principal
directions of the first-pass correction field:

``` r

bs <- fit2$two_pass$batch_subspace
cat("Batch subspace rank (pass 1):        ", bs$rank, "\n")
#> Batch subspace rank (pass 1):         2
cat("Variance explained by leading directions:", round(bs$variance_explained, 3), "\n")
#> Variance explained by leading directions: 0.931
```

This is a read-only diagnostic. Rank 1 indicates that pass-1 node shifts
are nearly collinear — which is expected when a single global offset
dominates. Higher rank means the first-pass correction has more
heterogeneous directions across the codebook, pointing toward a
spatially varied batch effect.

The subspace is **not** used for correction. Total shift = residual
pass-2 shift + global shift `g`. The diagnostic may conflate batch
effects with biology (first-pass OT sees the original unmodified
codebook, not controlled anchors), so it should be read descriptively,
not as a geometric estimate of the true batch direction. For an
anchor-based analogue that does serve as a geometric estimate, see
[`vignette("anchor-samples", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/anchor-samples.md).

The `variance_threshold` argument controls rank selection for this
diagnostic:

``` r

fit2b <- somalign_fit_two_pass(query, reference, variance_threshold = 1.0)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.41); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.08); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
cat("Rank at variance_threshold = 1.0:", fit2b$two_pass$batch_subspace$rank, "\n")
#> Rank at variance_threshold = 1.0: 4

# Correction is unchanged regardless of variance_threshold
cat("Node shifts identical:", isTRUE(all.equal(fit2$node_shifts, fit2b$node_shifts)), "\n")
#> Node shifts identical: TRUE
```

## Comparison with single-pass

For reference, a plain
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
at `epsilon_local` sees the full uncorrected offset in one pass:

``` r

fit1 <- somalign_fit(query, reference, epsilon = 0.1)
#> somalign_fit: 4 query node(s) have match_mass_ratio > 1 (max 1.10); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.

cat("Mean node shift -- two-pass: ",
    round(mean(sqrt(rowSums(fit2$node_shifts^2))), 4), "\n")
#> Mean node shift -- two-pass:  0.2982
cat("Mean node shift -- one-pass: ",
    round(mean(sqrt(rowSums(fit1$node_shifts^2))), 4), "\n")
#> Mean node shift -- one-pass:  0.209
```

Both return `somalign_fit` objects with an identical downstream
interface:

``` r

res2 <- somalign_results(fit2)
res1 <- somalign_results(fit1)
cat("Column names match:", identical(names(res2), names(res1)), "\n")
#> Column names match: TRUE
cat("Rows:", nrow(res2), "\n")
#> Rows: 120

diag2 <- somalign_diagnostics(fit2)
cat("Solver converged:", diag2$solver$converged, "\n")
#> Solver converged: TRUE
cat("Transport mass:  ", round(diag2$ot$transport_mass, 4), "\n")
#> Transport mass:   1.0527
```

## Label-guided alignment

Both
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
and
[`somalign_fit_two_pass()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_two_pass.md)
accept `label_guided = TRUE`, which applies a large cost penalty to node
pairs with discordant dominant labels. When reference labels are known
and the batch does not scramble population boundaries, this concentrates
OT transport onto same-label pairs and reduces spurious cross-cluster
mass.

`label_guided = TRUE` requires `query$label_prob` to be non-NULL. The
standard
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
workflow using
[`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html) leaves
this field as NULL; it is populated only when a
[`kohonen::supersom()`](https://rdrr.io/pkg/kohonen/man/supersom.html)
or [`kohonen::xyf()`](https://rdrr.io/pkg/kohonen/man/supersom.html) SOM
with a label code layer is passed as `som_query`. Build the query SOM
with a one-hot label layer, supply it via `som_query`, then set
`label_guided = TRUE`:

``` r

# One-hot label matrix for new_data (requires knowing new-batch labels)
Y <- model.matrix(~ old_labels - 1)
colnames(Y) <- levels(factor(old_labels))

# Train a supersom with the label layer
new_scaled <- scale(new_data, center = reference$center, scale = reference$scale)
query_ss   <- kohonen::supersom(
  list(new_scaled, Y),
  grid = kohonen::somgrid(3, 3, "hexagonal"),
  rlen = 30,
  keep.data = TRUE
)

query_labeled <- somalign_query(
  new_data, reference,
  som_query      = query_ss,
  codebook_space = "reference_scaled"
)

fit_lg <- somalign_fit_two_pass(
  query_labeled, reference,
  label_guided = TRUE
)
```

This works equally for
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
If labels are noisy or the batch shifts populations across label
boundaries, the cost penalty can destabilise the OT plan; check
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
and compare against
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
before using label-transferred assignments in downstream analyses.

## When to use `somalign_fit_two_pass()`

Two-pass alignment is most useful when the global batch offset is large
relative to `epsilon_local` — large enough that a single pass at low
epsilon struggles to route mass across the full displacement, while a
single pass at high epsilon loses population-level resolution. The
explicit global subtraction removes that tension. When the batch shift
is small or the OT problem is well-conditioned at a single epsilon,
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
is simpler and equally effective.

If remeasured control samples are available,
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
with `correction = "subspace"` offers geometrically grounded
signal-preserving correction as an alternative to the two-pass
decomposition.

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
