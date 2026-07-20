# Anchor-regularized alignment with remeasured samples

## When to use this

Reach for
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
when a set of quality-control samples was physically run in **both**
batches: control cell lines, reference bead populations, or
proficiency-panel specimens that appear in the old batch (which trained
the reference SOM) and the new batch (which trained the query SOM). To
use it you need two matched matrices with the same markers, `anchor_old`
and `anchor_new`, holding the same QC material measured under each batch
condition.

Those paired measurements pin down the batch shift directly. Wherever an
anchor sample lands in the reference, its query-batch measurement tells
you how far that region moved.
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
folds that information into the optimal-transport solve, so
anchor-supported node pairs are cheaper to route mass through. Nodes
with no anchor coverage keep the standard
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
cost, so anchors help where you have them and change nothing where you
do not.

Anchors shape the *correction* (the per-node shift), which `somalign`
treats as a diagnostic rather than a product. They still help label
transfer indirectly, because the same anchor cost bonus biases the
transport plan that labels are read from. The corrected coordinates
themselves can over-merge populations (see
[`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md))
and are not a batch-corrected expression matrix.

### How the anchor pairs become transport costs

For each anchor pair, the old-batch measurement is projected onto the
query SOM and the new-batch measurement onto the reference SOM. The
query SOM was trained on new-batch data, so projecting the old-batch
anchor onto it identifies which query node that anchor occupied before
the shift. Projecting the new-batch anchor onto the reference SOM
identifies the matching reference node after the shift. Counting these
pairs gives a node-pair matrix: how many anchors link query node *k* to
reference node *l*. Pairs with anchor support get reduced transport
cost, so the plan preferentially routes mass through them.

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

The new batch carries a batch effect that differs by marker, realistic
for a reagent lot change or inter-instrument drift where each channel
shifts by its own amount. A small per-cell jitter reflects run-to-run
variation, so the anchor displacements span a low-rank subspace rather
than one identical vector (this matters for the subspace correction
below).

``` r

batch_shift <- c(0.6, -0.4, 0.5, 0.0, 0.3)          # per-marker technical offset
new_data <- old_data +
  matrix(batch_shift, nrow = nrow(old_data), ncol = p, byrow = TRUE) +
  matrix(rnorm(nrow(old_data) * p, sd = 0.15), ncol = p)   # per-cell jitter

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

Suppose 25 samples from a QC pool (the same biological material measured
under both batch conditions) are spread across both populations.
`anchor_old` and `anchor_new` are the same rows measured in each batch.

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
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.30); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
fit_anc   <- somalign_fit_anchored(query, reference,
                                    anchor_old = anchor_old,
                                    anchor_new = anchor_new,
                                    rho_anchor = 1.5)
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
```

Printing either object gives a concise summary:

``` r

fit_plain
#> <somalign_fit>
#>   label transfer: disabled (reference has no labels)
#>   solver: internal  |  query nodes: 16  |  reference nodes: 16  |  transport mass: 1.144
fit_anc
#> <somalign_anchored_fit>
#>   label transfer: disabled (reference has no labels)
#>   solver: internal  |  query nodes: 16  |  reference nodes: 16  |  transport mass: 1.159
#>   anchors: 25 (56.2% node coverage) -- correction is a diagnostic, not a corrected-expression product
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
#> $correction
#> [1] "cost_bonus"
#> 
#> $nodes_covered
#> [1] 9
#> 
#> $coverage_fraction
#> [1] 0.5625
#> 
#> $batch_subspace
#> NULL
#> 
#> $variance_threshold
#> [1] 0.9
#> 
#> $displacements
#>               CD1        CD2         CD3           CD4           CD5
#>  [1,] -0.32948579 0.24747754 -0.20811356 -0.0787562200 -0.0823889701
#>  [2,] -0.29779831 0.26268504 -0.09190164  0.2140135925 -0.1661392463
#>  [3,] -0.21514604 0.12578006 -0.26260956 -0.1464489718 -0.0749789151
#>  [4,] -0.36031636 0.27339547 -0.13628896  0.0780024293 -0.0001725585
#>  [5,] -0.35110492 0.09941999 -0.15139468  0.0874705423  0.0221494496
#>  [6,] -0.23812394 0.25251676 -0.27277401  0.0244089884 -0.2069312954
#>  [7,] -0.48112593 0.21813122 -0.28631511 -0.0469622657 -0.1498354072
#>  [8,] -0.06199709 0.25376616 -0.20796482  0.0607195338 -0.1721443307
#>  [9,] -0.33453366 0.11542366 -0.11898160  0.1337829235  0.0043469955
#> [10,] -0.25341484 0.28997653 -0.06895253  0.0349333936 -0.1338518379
#> [11,] -0.30190450 0.13986891 -0.13605760 -0.0234081872 -0.1514135114
#> [12,] -0.21246930 0.12268142 -0.18556322 -0.0616977375 -0.2635092129
#> [13,] -0.27337707 0.18531817 -0.25966850  0.0706363342 -0.2054948563
#> [14,] -0.18080448 0.21788065 -0.12164558  0.0745978634 -0.0516513986
#> [15,] -0.22999057 0.25944304 -0.21189366 -0.0969612572 -0.0756928273
#> [16,] -0.27626718 0.14959996 -0.08247805  0.0824329907 -0.2173444277
#> [17,] -0.32588432 0.23153357 -0.23594851 -0.0186016396 -0.1366205239
#> [18,] -0.38900213 0.13184508 -0.16473659  0.1111233033 -0.2584214396
#> [19,] -0.25732940 0.22320138 -0.13327184  0.0005092108 -0.3398426151
#> [20,] -0.27271323 0.20387518 -0.24651942 -0.0404862614 -0.1586898798
#> [21,] -0.27140128 0.25266406 -0.36303814 -0.0792324542 -0.1992333866
#> [22,] -0.25258361 0.21721280 -0.12308932  0.0798243950 -0.1525217546
#> [23,] -0.26588694 0.16612937 -0.17310313  0.0531524997 -0.1037191411
#> [24,] -0.33938973 0.33836623 -0.22082513 -0.0638744551 -0.0925247562
#> [25,] -0.41450561 0.26612702 -0.18486929 -0.0420923984 -0.2220232949
#> 
#> $feature_weights
#> NULL
```

`nodes_covered` counts distinct query SOM nodes onto which at least one
anchor’s old-batch measurement was projected. `coverage_fraction` is
that count over the total number of query nodes. Nodes with no anchor
coverage are solved with unchanged costs, so they behave exactly as they
would under
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## Effect on node correction vectors

Each correction vector points from a query SOM node toward its
barycentric target in the reference, so a node whose OT correspondence
changed shows a different norm. Comparing per-node correction magnitudes
is a quick way to see the anchor effect:

``` r

plain_norms <- sqrt(rowSums(fit_plain$node_shifts^2))
anc_norms   <- sqrt(rowSums(fit_anc$node_shifts^2))

cat("Mean correction norm  plain:    ", round(mean(plain_norms), 4), "\n")
#> Mean correction norm  plain:     0.4393
cat("Mean correction norm  anchored: ", round(mean(anc_norms),   4), "\n")
#> Mean correction norm  anchored:  0.4303
cat("Max  correction norm  plain:    ", round(max(plain_norms),  4), "\n")
#> Max  correction norm  plain:     0.6898
cat("Max  correction norm  anchored: ", round(max(anc_norms),    4), "\n")
#> Max  correction norm  anchored:  0.6906
```

When anchors consistently confirm the same direction of displacement,
the anchored fit produces more coherent correction vectors. The plan
concentrates mass on anchor-supported routes, which reduces the
influence of diffuse entropic transport on uncovered node pairs.

## Downstream results

Both fit objects expose the same downstream interface.
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
returns the standard per-sample projection data frame:

``` r

res_plain <- somalign_results(fit_plain)
res_anc   <- somalign_results(fit_anc)

# Fraction of samples landing outside reference distance thresholds
cat("Outside fraction (plain):    ",
    mean(res_plain$final_status == "outside_reference"), "\n")
#> Outside fraction (plain):     0.4333333
cat("Outside fraction (anchored): ",
    mean(res_anc$final_status   == "outside_reference"), "\n")
#> Outside fraction (anchored):  0.4333333
```

[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
works identically on anchored fits:

``` r

diag <- somalign_diagnostics(fit_anc)
cat("Solver converged:", diag$solver$converged, "\n")
#> Solver converged: TRUE
cat("Transport mass:  ", round(diag$ot$transport_mass, 4), "\n")
#> Transport mass:   1.1589
```

## Tuning `rho_anchor`

`rho_anchor` controls how much the anchor bonus reduces the normalized
OT cost. At `rho_anchor = 0` the bonus vanishes and the result equals
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
As `rho_anchor` grows, anchor-supported pairs get progressively cheaper.
Once the bonus exceeds the normalized cost for a pair, that pair’s
effective cost clamps to zero; the plan then spreads mass among all
zero-cost pairs by entropic regularization rather than by differential
anchor support, so very large values stop concentrating mass further.

A practical upper bound is `rho_anchor * max_A / n_anchors <= 1`, where
`max_A` is the highest anchor count for any single node pair. Values
from 0.5 to 2 are a reasonable starting range. When coverage is sparse
(below 30% of nodes) or anchor pairs are unevenly distributed, a smaller
value avoids over-weighting a few pairs.

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
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
do.call(rbind, sweep_results)
#>   rho_anchor transport_mass mean_corr_norm
#> 1        0.5        1.15210        0.43308
#> 2        1.0        1.15589        0.43114
#> 3        1.5        1.15894        0.43026
#> 4        2.0        1.16115        0.43055
#> 5        4.0        1.16820        0.43549
```

For a more targeted diagnostic, compare the
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
node-level `match_fraction` and `correction_norm` across values, and use
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
to check whether corrected node assignments move in a direction
consistent with your biological expectations.

## Signal-preserving correction

The default `correction = "cost_bonus"` applies the OT-derived shift
across all markers. For most reagent- or calibration-driven batch
effects this is correct: the displacement is genuinely technical and
should be removed. The assumption breaks when the new batch contains
populations or activation states absent from the reference. A full-space
shift then pushes those cells toward reference coordinates and erases
their distinguishing biology.

Each anchor pair links the same biological unit measured in both
batches, so the row vector `anchor_old - anchor_new` observes the batch
displacement with no biological contribution. An SVD of the `n_anchors`
by `p` displacement matrix isolates the batch directions; call this
`V_batch`. Restricting node shifts to `V_batch` corrects only the batch
component and leaves biology orthogonal to it untouched.

``` r

fit_sub <- somalign_fit_anchored(
  query, reference,
  anchor_old = anchor_old,
  anchor_new = anchor_new,
  correction = "subspace"
)
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.30); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
cat("Correction mode:  ", fit_sub$anchors$correction, "\n")
#> Correction mode:   subspace
cat("Batch subspace V: ", nrow(fit_sub$anchors$batch_subspace$V), "rows x",
                          ncol(fit_sub$anchors$batch_subspace$V), "cols\n")
#> Batch subspace V:  5 rows x 2 cols
cat("Rank:             ", fit_sub$anchors$batch_subspace$rank, "\n")
#> Rank:              2
cat("Variance explained:", round(fit_sub$anchors$batch_subspace$variance_explained, 3), "\n")
#> Variance explained: 0.932
```

The `batch_subspace` element carries `V` (p by rank), `rank`, and
`variance_explained`. The default `variance_threshold = 0.9` selects the
smallest rank r whose top r squared singular values account for at least
90% of the total. A pure uniform offset collapses to rank 1; a
heterogeneous batch effect, like the per-marker shift plus jitter here,
occupies more dimensions. Raising the threshold to 1.0 keeps every
direction, so the rank grows without changing the anchor data:

``` r

fit_sub_full <- somalign_fit_anchored(
  query, reference,
  anchor_old        = anchor_old,
  anchor_new        = anchor_new,
  correction        = "subspace",
  variance_threshold = 1.0
)
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.30); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
cat("Rank at threshold = 0.9:", fit_sub$anchors$batch_subspace$rank, "\n")
#> Rank at threshold = 0.9: 2
cat("Rank at threshold = 1.0:", fit_sub_full$anchors$batch_subspace$rank, "\n")
#> Rank at threshold = 1.0: 5
```

`correction = "both"` applies the anchor cost bonus to the OT solve and
restricts the resulting shifts to `V_batch`:

``` r

fit_both <- somalign_fit_anchored(
  query, reference,
  anchor_old = anchor_old,
  anchor_new = anchor_new,
  rho_anchor = 1.5,
  correction = "both"
)
#> somalign_fit: 13 query node(s) have match_mass_ratio > 1 (max 1.29); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
cat("Correction mode:", fit_both$anchors$correction, "\n")
#> Correction mode: both
cat("Coverage:       ", round(fit_both$anchors$coverage_fraction, 3), "\n")
#> Coverage:        0.562
cat("Rank:           ", fit_both$anchors$batch_subspace$rank, "\n")
#> Rank:            2
```

`rho_anchor` has no effect under `correction = "subspace"`, since there
is no cost bonus to apply; under `"both"` the same tuning guidance as
`"cost_bonus"` holds. For a simple uniform offset all three modes give
similar corrected positions. The difference matters when anchor coverage
is uneven or the query carries populations the reference does not
represent.

## Corrected marker expression for downstream analysis

Label transfer is the primary product, but some downstream steps, such
as a UMAP or a differential-expression contrast, want a corrected
expression matrix rather than a per-cell label.
[`somalign_correct_expression()`](https://mdmanurung.github.io/somalign/reference/somalign_correct_expression.md)
returns one: a cells-by-markers matrix, smoothed across each cell’s
nearest SOM nodes. For an anchored fit with `correction = "subspace"` or
`"both"` (like `fit_sub` here), the correction is confined to the
anchor-estimated batch subspace, so variation orthogonal to the batch
direction stays intact. A two-pass fit carries no anchor subspace;
[`somalign_correct_expression()`](https://mdmanurung.github.io/somalign/reference/somalign_correct_expression.md)
then applies its full two-pass shift instead. Either kind of fit works;
a plain
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
carries no batch correction and is rejected.

``` r

expr_corrected <- somalign_correct_expression(fit_sub)   # raw marker units
head(expr_corrected[, 1:4])
```

Treat this output as an aid for visualisation and differential
expression, not the primary product. The correction preserves orthogonal
structure, but within the batch subspace it reduces rather than fully
removes the distance between populations. Run
`somalign_topology_audit(fit_sub)` first to confirm the correction is
warranted, because populations sitting close together along the batch
direction can be drawn toward one another. To compare cell-type
composition or abundance across batches, prefer the direct projection
columns from
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
with a compositional (centred log-ratio) transform.

## Choosing a mode, and when to skip anchors

[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
fits when a set of samples (QC pools, reference standards, proficiency
panels) ran in both batches, giving ground-truth node-pair
correspondences. Anchor pairs need not cover the full codebook; even
sparse coverage constrains the plan for the covered pairs, while
uncovered nodes fall back to the standard objective.

If no remeasured samples exist, use
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
directly. Fabricated or imputed anchor pairs add no information and will
misguide the transport. When query-only biology must be preserved, use
`correction = "subspace"` or `"both"`; when a full-space shift is
appropriate (the batch dominates and no signal-preservation is needed),
the default `"cost_bonus"` is enough. See
[`vignette("two-pass", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/two-pass.md)
for the two-pass decomposition, which handles large global offsets
without remeasured anchors.

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
    #> [1] somalign_0.99.5  kohonen_3.0.13   BiocStyle_2.40.0
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
