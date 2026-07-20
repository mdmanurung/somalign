# API smoke test: all public functions

This vignette exercises every public function in the package so that
build-time vignette rendering catches breaking changes or silent
regressions. It is not a tutorial (see the other vignettes for annotated
workflows). Each section calls the function, checks the return type, and
prints a key field to confirm the object was constructed correctly.

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
#>   label transfer: 100.0% of cells accepted across 2 class(es); median confidence 1.00, median margin 1.00
#>   solver: internal  |  query nodes: 9  |  reference nodes: 9  |  transport mass: 1.104
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
#> Rows: 60   Columns: 23

res_aug <- somalign_results(fit, data = new_data)
stopifnot(ncol(res_aug) > ncol(res))
```

## `somalign_soft_labels`

Computes per-cell soft label probabilities by smoothing over each cell’s
nearest reference SOM nodes.

``` r

soft <- somalign_soft_labels(fit, k = 3)
stopifnot(inherits(soft, "somalign_soft_labels"))
stopifnot(nrow(soft) == nrow(new_data))
cat("Soft label columns:", paste(colnames(soft), collapse = ", "), "\n")
#> Soft label columns: pop_A, pop_B
```

## `somalign_soft_frequencies`

Aggregates soft labels by sample or another grouping variable.

``` r

sample_group <- rep(c("sample_1", "sample_2"), length.out = nrow(new_data))
freq <- somalign_soft_frequencies(fit, group = sample_group, k = 3)
stopifnot(inherits(freq, "somalign_soft_frequencies"))
stopifnot(nrow(freq) == 2L)
cat("Frequency row sums:", paste(round(rowSums(freq), 3), collapse = ", "), "\n")
#> Frequency row sums: 1, 1
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

## `somalign_epsilon_sweep`

Runs an OT-only sweep over epsilon values and reports plan-geometry
diagnostics without re-projecting cells.

``` r

eps <- somalign_epsilon_sweep(
  qry, ref,
  epsilon_grid = c(0.05, 0.1, 0.2),
  solver = "internal"
)
stopifnot(inherits(eps, "somalign_epsilon_sweep"))
stopifnot(nrow(eps$table) == 3L)
cat("epsilon_rec:", round(eps$epsilon_rec, 4), "\n")
#> epsilon_rec: 0.03
```

## `somalign_select_epsilon`

Selects an epsilon from the sweep curve using an unsupervised criterion.

``` r

sel <- somalign_select_epsilon(
  qry, ref,
  epsilon = c(0.05, 0.1, 0.2),
  solver = "internal"
)
stopifnot(inherits(sel, "somalign_epsilon_selection"))
cat("selected epsilon:", round(sel$selected_epsilon, 4), "\n")
#> selected epsilon: 0.03
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
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.27); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.22); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
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

## `somalign_topology_audit`

Checks whether the barycentric correction collapses SOM topology.

``` r

topo <- somalign_topology_audit(fit, nodes = "all")
#> Warning: topology_warning: corrected codebook has 4 H0 component(s) vs 5 in
#> query (delta = -1; populations may have been merged/erased). Inspect
#> fit$diagnostics$topology for details.
stopifnot(inherits(topo, "somalign_topology"))
topo
#> <somalign_topology>
#>   threshold: 0.4515 (auto)
#>   H0 components  query: 5  corrected: 4  reference: 3
#>   topology_delta: -1   warning: TRUE
```

## `somalign_worst_nodes`

Returns the query nodes with the lowest match fraction.

``` r

worst <- somalign_worst_nodes(fit, n = 3)
stopifnot(is.data.frame(worst))
stopifnot(nrow(worst) == 3L)
worst[, c("query_node", "match_fraction")]
#>    query_node match_fraction
#> V1          1              0
#> V5          5              0
#> V8          8              0
```

## `somalign_check_codebook_alignment`

Preflight check for query/reference codebook coordinate compatibility.

``` r

chk <- somalign_check_codebook_alignment(
  qry$codebook, ref,
  query_masses = qry$node_masses,
  stop_if_critical = FALSE
)
stopifnot(inherits(chk, "somalign_codebook_check"))
cat("verdict:", chk$verdict, "\n")
#> verdict: pass
```

## `somalign_label_metrics`

Scores predicted labels against ground truth: accuracy, macro-F1,
multiclass MCC, per-class statistics, and the confusion matrix.

``` r

truth <- rep(c("A", "B"), each = 10)
pred  <- truth; pred[1] <- "B"
metrics <- somalign_label_metrics(pred, truth)
metrics
#> <somalign_label_metrics>
#>   accuracy = 0.9500  macro_f1 = 0.9499  MCC = 0.9045
#>   scored = 20  coverage = 100.0%  accuracy_all = 0.9500
```

## `somalign_calibration`

Bins predictions by a confidence score and compares mean confidence to
empirical accuracy, returning the expected calibration error, Brier
score, and the coverage (fraction of non-abstained predictions actually
scored).

``` r

calib <- somalign_calibration(runif(50), runif(50) < 0.5)
calib
#> <somalign_calibration>
#>   ECE = 0.2827  MCE = 0.7356  Brier = 0.3287  (scored n = 50, coverage = 100.0%)
#>   reliability (score_mean -> accuracy, n):
#>     0.05 -> 0.57  (7)
#>     0.16 -> 0.67  (3)
#>     0.26 -> 1.00  (1)
#>     0.35 -> 0.43  (7)
#>     0.44 -> 0.25  (8)
#>     0.63 -> 0.83  (6)
#>     0.77 -> 0.67  (3)
#>     0.87 -> 0.67  (6)
#>     0.94 -> 0.56  (9)
```

## `somalign_cross_validate`

Held-out k-fold cross-validation of label transfer.

``` r

cv <- somalign_cross_validate(
  old_data, labels_old, grid = grid_small, k = 2, rlen = 15
)
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.30); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
stopifnot(inherits(cv, "somalign_cross_validation"))
cv
#> <somalign_cross_validation> (2 folds)
#>   pooled: accuracy = 1.0000  macro_f1 = 1.0000  MCC = 1.0000  coverage = 100.0%
#>   calibration: ECE = 0.0000  Brier = 0.0000  (scored on 100.0% of predictions)
```

## `somalign_tune`

Selects transport-plan parameters by cross-validated label accuracy.

``` r

tuned <- somalign_tune(
  old_data, labels_old, grid = grid_small,
  param_grid = data.frame(epsilon = c(0.1, 0.2)), k = 2, rlen = 15
)
stopifnot(inherits(tuned, "somalign_tune"))
tuned$grid
#>   epsilon rho_query rho_ref diagonal_boost feature_weights accuracy macro_f1
#> 1     0.1         1       1              0            none        1        1
#> 2     0.2         1       1              0            none        1        1
#>   mcc coverage          ece
#> 1   1        1 0.000000e+00
#> 2   1        1 1.835643e-12
```

## `somalign_anchor_benefit`

Quantifies the label-transfer lift from anchor samples across a
`rho_anchor` grid.

``` r

anc_i <- c(1:10, 31:40)
ab <- somalign_anchor_benefit(
  qry, ref, labels_old[seq_len(nrow(new_data))],
  anchor_old = old_data[anc_i, ], anchor_new = new_data[anc_i, ],
  rho_grid = c(0, 1)
)
stopifnot(inherits(ab, "somalign_anchor_benefit"))
ab$grid
#>   rho_anchor accuracy macro_f1 mcc coverage          ece
#> 1          0        1        1   1        1 1.998401e-15
#> 2          1        1        1   1        1 1.443290e-15
```

## `somalign_correct_expression`

Returns an auxiliary cell-by-marker corrected expression matrix for
anchored subspace or two-pass fits.

``` r

expr <- somalign_correct_expression(fit_sub, k = 3)
stopifnot(inherits(expr, "somalign_corrected_expression"))
stopifnot(identical(dim(expr), dim(new_data)))
cat("Corrected expression dimensions:", paste(dim(expr), collapse = " x "), "\n")
#> Corrected expression dimensions: 60 x 4
```

## `somalign_subspace_sensitivity`

Bootstraps the anchor-derived batch subspace and reports correction
uncertainty.

``` r

sens <- somalign_subspace_sensitivity(fit_sub, n_boot = 5L)
stopifnot(inherits(sens, "somalign_subspace_sensitivity"))
cat("Subspace rank:", sens$subspace_rank, "\n")
#> Subspace rank: 1
```

## `somalign_exclusion_test`

Permutation diagnostic for structured residual signal outside the anchor
batch subspace.

``` r

excl <- somalign_exclusion_test(fit_sub, n_perm = 9L)
stopifnot(inherits(excl, "somalign_exclusion_test"))
cat("Exclusion verdict:", excl$verdict, "\n")
#> Exclusion verdict: pass
```

## Plot functions

All plotting helpers return `ggplot` objects.

``` r

plot_list <- list(
  somalign_plot_mass_balance(fit),
  somalign_plot_match_fraction(fit),
  somalign_plot_correction(fit),
  somalign_plot_outside_fraction(fit),
  somalign_plot_label_confusion(fit),
  somalign_plot_codebook_range(chk),
  somalign_plot_marker_distributions(qry, reference = ref)
)
stopifnot(all(vapply(plot_list, inherits, logical(1), what = "ggplot")))
```

The epsilon sweep object also has a plot method:

``` r

eps_plot <- plot(eps)
stopifnot(inherits(eps_plot, "ggplot"))
```

## Print and summary methods

The main classes have print methods that return the object invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) on a fit reports the
label-transfer breakdown.

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
#>   label transfer: 100.0% of cells accepted across 2 class(es); median confidence 1.00, median margin 1.00
#>   solver: internal  |  query nodes: 9  |  reference nodes: 9  |  transport mass: 1.104
invisible(print(fit_anc))
#> <somalign_anchored_fit>
#>   label transfer: 100.0% of cells accepted across 2 class(es); median confidence 1.00, median margin 1.00
#>   solver: internal  |  query nodes: 9  |  reference nodes: 9  |  transport mass: 1.123
#>   anchors: 20 (66.7% node coverage) -- correction is a diagnostic, not a corrected-expression product
invisible(print(eps))
#> somalign epsilon sweep [3 points]
#>   epsilon range      : 0.05 - 0.2
#>   critical epsilon   : 0.1
#>   recommended epsilon (0.3x critical): 0.03
#>   cost_scale         : 3.755
invisible(print(metrics))
#> <somalign_label_metrics>
#>   accuracy = 0.9500  macro_f1 = 0.9499  MCC = 0.9045
#>   scored = 20  coverage = 100.0%  accuracy_all = 0.9500
invisible(print(calib))
#> <somalign_calibration>
#>   ECE = 0.2827  MCE = 0.7356  Brier = 0.3287  (scored n = 50, coverage = 100.0%)
#>   reliability (score_mean -> accuracy, n):
#>     0.05 -> 0.57  (7)
#>     0.16 -> 0.67  (3)
#>     0.26 -> 1.00  (1)
#>     0.35 -> 0.43  (7)
#>     0.44 -> 0.25  (8)
#>     0.63 -> 0.83  (6)
#>     0.77 -> 0.67  (3)
#>     0.87 -> 0.67  (6)
#>     0.94 -> 0.56  (9)
invisible(print(cv))
#> <somalign_cross_validation> (2 folds)
#>   pooled: accuracy = 1.0000  macro_f1 = 1.0000  MCC = 1.0000  coverage = 100.0%
#>   calibration: ECE = 0.0000  Brier = 0.0000  (scored on 100.0% of predictions)
invisible(print(tuned))
#> <somalign_tune> objective = mcc over 2 combination(s)
#>   best:
#>     epsilon=0.1 rho_query=1 rho_ref=1 diagonal_boost=0 fw=none
#>     accuracy=1.0000 macro_f1=1.0000 MCC=1.0000 coverage=100.0% ECE=0.0000
invisible(print(ab))
#> <somalign_anchor_benefit> objective = mcc over 2 rho value(s)
#>   baseline (rho=0): mcc = 1.0000  coverage = 100.0%
#>   best: rho=0  mcc = 1.0000  (lift +0.0000)  coverage = 100.0%  ECE = 0.0000
invisible(print(topo))
#> <somalign_topology>
#>   threshold: 0.4515 (auto)
#>   H0 components  query: 5  corrected: 4  reference: 3
#>   topology_delta: -1   warning: TRUE
invisible(print(sens))
#> <somalign_subspace_sensitivity>
#>   n_anchors = 20   rank = 1   n_boot = 5   conf_level = 0.95
#>   median principal angle: 0.0 deg
#>   median tipping angle:    48.0 deg (n=6 allowed nodes)
invisible(print(excl))
#> <somalign_exclusion_test>
#>   sv_observed = 0.0000   p = 0.4444   verdict = pass
#>   null sv (2.5%/50%/97.5%): 0.0000 / 0.0000 / 0.0000
#>   rank = 1   n_anchors = 20   n_features = 4
invisible(summary(fit))
#> <somalign_fit> label-transfer summary
#>   cells: 60  |  accepted: 60 (100.0%)  |  classes: 2
#>   confidence quartiles (accepted): 1.00 / 1.00 / 1.00
#>   median margin (accepted): 1.00
#>   accepted class distribution:
#>     pop_A                30
#>     pop_B                30
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
    #> [1] somalign_0.99.5  kohonen_3.0.13   BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.2.1        
    #>  [4] compiler_4.6.1      BiocManager_1.30.27 tidyselect_1.2.1   
    #>  [7] Rcpp_1.1.2          jquerylib_0.1.4     scales_1.4.0       
    #> [10] yaml_2.3.12         fastmap_1.2.0       ggplot2_4.0.3      
    #> [13] R6_2.6.1            generics_0.1.4      knitr_1.51         
    #> [16] htmlwidgets_1.6.4   tibble_3.3.1        bookdown_0.47      
    #> [19] desc_1.4.3          bslib_0.11.0        pillar_1.11.1      
    #> [22] RColorBrewer_1.1-3  rlang_1.3.0         cachem_1.1.0       
    #> [25] xfun_0.60           fs_2.1.0            sass_0.4.10        
    #> [28] S7_0.2.2            otel_0.2.0          viridisLite_0.4.3  
    #> [31] cli_3.6.6           pkgdown_2.2.1       withr_3.0.3        
    #> [34] magrittr_2.0.5      digest_0.6.39       grid_4.6.1         
    #> [37] lifecycle_1.0.5     vctrs_0.7.3         evaluate_1.0.5     
    #> [40] glue_1.8.1          farver_2.1.2        rmarkdown_2.31     
    #> [43] pkgconfig_2.0.3     tools_4.6.1         htmltools_0.5.9
