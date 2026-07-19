# Validating and tuning label transfer

The product of a `somalign` fit is a transferred label and a confidence
for every query cell. This vignette covers how to check that those
labels are correct and that the confidence is honest, and how to choose
transport parameters that maximise label accuracy rather than an
unsupervised proxy. The corrected coordinates are a separate, diagnostic
output and are not used here; see the algorithm vignette for why they
can over-merge populations.

## Scoring labels against a known truth

[`somalign_label_metrics()`](https://mdmanurung.github.io/somalign/reference/somalign_label_metrics.md)
scores predicted labels against ground truth. It reports overall
accuracy, macro-averaged F1, the multiclass Matthews correlation
coefficient (Gorodkin’s ), per-class precision, recall and F1, and the
confusion matrix. Metrics are computed on the accepted predictions;
`coverage` records what fraction of cells that is.

``` r

library(somalign)

truth <- rep(c("Tcell", "Bcell", "Mono"), each = 20)
pred  <- truth
pred[c(1, 25, 45)] <- c("Bcell", "Tcell", "Tcell")   # three errors of 60

m <- somalign_label_metrics(pred, truth)
m
#> <somalign_label_metrics>
#>   accuracy = 0.9500  macro_f1 = 0.9504  MCC = 0.9254
#>   scored = 60  coverage = 100.0%  accuracy_all = 0.9500
m$per_class
#>   class precision recall        f1 support
#> 1 Bcell 0.9500000   0.95 0.9500000      20
#> 2  Mono 1.0000000   0.95 0.9743590      20
#> 3 Tcell 0.9047619   0.95 0.9268293      20
```

The `accepted` argument marks which predictions to score; cells that are
not accepted count as abstentions and are reported through `coverage`
and `accuracy_all` rather than as errors.

## Is the confidence honest?

A confidence of 0.8 should mean the label is right about 80% of the
time.
[`somalign_calibration()`](https://mdmanurung.github.io/somalign/reference/somalign_calibration.md)
bins predictions by their confidence and compares the mean confidence in
each bin to the empirical accuracy there, summarising the gap as the
expected calibration error (ECE) and the Brier score.

``` r

set.seed(1)
score   <- runif(500)
correct <- runif(500) < score          # calibrated by construction
somalign_calibration(score, correct)
#> <somalign_calibration>
#>   ECE = 0.0509  MCE = 0.1869  Brier = 0.1729  (n = 500)
#>   reliability (score_mean -> accuracy, n):
#>     0.06 -> 0.05  (38)
#>     0.15 -> 0.16  (61)
#>     0.25 -> 0.22  (50)
#>     0.35 -> 0.32  (57)
#>     0.45 -> 0.52  (64)
#>     0.54 -> 0.52  (42)
#>     0.65 -> 0.47  (43)
#>     0.75 -> 0.82  (50)
#>     0.85 -> 0.79  (43)
#>     0.95 -> 0.98  (52)
```

A well-calibrated model has mean confidence close to accuracy in every
bin, so its ECE is near zero. A model that reports 0.99 on cells it gets
right only half the time has an ECE near 0.49.

## Held-out cross-validation

The functions above need ground-truth labels, which a real query batch
does not carry.
[`somalign_cross_validate()`](https://mdmanurung.github.io/somalign/reference/somalign_cross_validate.md)
supplies them by construction: it splits the labelled reference into
folds, trains a reference on each training split, projects the held-out
split, and transfers labels back. Every held-out cell has a real label,
so the pooled score estimates how well transfer generalises to cells the
reference never saw.

``` r

set.seed(1)
x <- rbind(
  matrix(rnorm(300 * 3, -3, 0.5), ncol = 3),
  matrix(rnorm(300 * 3,  3, 0.5), ncol = 3)
)
colnames(x) <- paste0("m", seq_len(3))
lab <- rep(c("low", "high"), each = 300)

cv <- somalign_cross_validate(
  x, lab, grid = kohonen::somgrid(3, 3, "hexagonal"), k = 3, rlen = 20
)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.16); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.17); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_fit: 6 query node(s) have match_mass_ratio > 1 (max 1.16); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
cv
#> <somalign_cross_validation> (3 folds)
#>   pooled: accuracy = 1.0000  macro_f1 = 1.0000  MCC = 1.0000  coverage = 100.0%
#>   calibration: ECE = 0.0000  Brier = 0.0000
cv$per_fold
#>   fold accuracy macro_f1 mcc coverage
#> 1    1        1        1   1        1
#> 2    2        1        1   1        1
#> 3    3        1        1   1        1
```

## Tuning the transport plan for accuracy

[`somalign_select_epsilon()`](https://mdmanurung.github.io/somalign/reference/somalign_select_epsilon.md)
picks a regularisation strength from the plan geometry alone, without
reference to labels. When labelled data is available,
[`somalign_tune()`](https://mdmanurung.github.io/somalign/reference/somalign_tune.md)
selects transport parameters by cross-validated label accuracy instead.
It sweeps a parameter grid, runs the same held-out cross-validation at
each point, and reports every metric so the operating point can be
chosen deliberately.

``` r

tuned <- somalign_tune(
  x, lab, grid = kohonen::somgrid(3, 3, "hexagonal"),
  param_grid = data.frame(epsilon = c(0.05, 0.1, 0.2)),
  k = 3, rlen = 20, metric = "mcc"
)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
tuned$grid
#>   epsilon rho_query rho_ref diagonal_boost feature_weights accuracy macro_f1
#> 1    0.05         1       1              0            none        1        1
#> 2    0.10         1       1              0            none        1        1
#> 3    0.20         1       1              0            none        1        1
#>   mcc coverage          ece
#> 1   1        1 0.000000e+00
#> 2   1        1 0.000000e+00
#> 3   1        1 4.850388e-10
```

The default objective is the multiclass MCC. Read it alongside
`coverage`: MCC and accuracy are scored on accepted cells only, so a
large `epsilon` can raise them by abstaining on hard cells. `macro_f1`
falls when rare classes are abstained away, which makes it a more
coverage-robust target for imbalanced data.

## When do anchor samples help?

Anchors are quality-control samples run in both batches. Through
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
they bias the transport plan toward anchor-supported node pairs.
[`somalign_anchor_benefit()`](https://mdmanurung.github.io/somalign/reference/somalign_anchor_benefit.md)
quantifies how much they improve label transfer by sweeping the anchor
strength `rho_anchor` against a known query label; `rho_anchor = 0` is
the no-anchor baseline.

``` r

set.seed(1)
shift <- 2
ref_x <- rbind(matrix(rnorm(120 * 3, -3, 0.5), ncol = 3),
               matrix(rnorm(120 * 3,  3, 0.5), ncol = 3))
colnames(ref_x) <- paste0("m", seq_len(3))
ref_lab <- rep(c("low", "high"), each = 120)
ref <- somalign_train_reference(ref_x, labels = ref_lab,
                                grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 20)

qry_x <- ref_x + shift
qry <- somalign_query(qry_x, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 20)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.

anc_idx <- c(1:30, 121:150)
ab <- somalign_anchor_benefit(
  qry, ref, ref_lab,
  anchor_old = ref_x[anc_idx, ],
  anchor_new = ref_x[anc_idx, ] + shift,
  rho_grid = c(0, 1, 10, 100)
)
ab$grid
#>   rho_anchor accuracy macro_f1 mcc coverage          ece
#> 1          0        1        1   1        1 1.432153e-04
#> 2          1        1        1   1        1 6.966324e-05
#> 3         10        1        1   1        1 2.906399e-06
#> 4        100        1        1   1        1 2.164835e-06
```

Anchors help only when the batch effect is population-specific and
severe enough that the plain plan mis-maps lineages, and then only at
`rho_anchor` well above 1. For mild or global batch shifts the plain
plan already succeeds, so the anchor lift is near zero.

## Quantifying abundance across batches

Label transfer gives each cell a cluster; comparing cluster abundance
across batches is a separate, compositional problem. Cluster frequencies
sum to one, so correlating raw proportions is dominated by the few
largest clusters and understates how well rarer clusters are recovered.
Quantify abundance with a centred log-ratio (CLR) transform of
per-sample cluster counts instead. The `crumblr` package computes the
CLR together with count-based precision weights suitable for
differential-abundance models:

``` r

library(crumblr)
counts <- table(sample_id, transferred_label)   # samples x clusters, integer
cobj <- crumblr(as.matrix(counts))               # CLR in cobj$E, weights in cobj$weights
```

Use the direct projection (`old_som_unit`) for composition, not the
corrected one; the CLR abundance profile is reproducible across batches,
whereas absolute per-cluster proportions are only approximate.

For the abundance step specifically, prefer **soft** projection over
hard nearest-node counts. Hard assignment sends each cell to a single
node and counts it once for that node’s cluster, so a cell near a
cluster boundary contributes a whole unit to one side and a small batch
shift can flip it — inflating the sampling variance of per-sample
cluster proportions.
[`somalign_soft_frequencies()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_frequencies.md)
instead spreads each cell over its nearest reference nodes, giving a
smoother per-sample frequency profile that reproduces better across
batches. Pass a node-to-cluster map through `node_groups` when the
grouping of interest is coarser than the reference labels:

``` r

# node2cluster: length-n_nodes map from reference node to metacluster.
# Use normalize = FALSE for crumblr, which models counts, not frequencies.
counts <- somalign_soft_frequencies(fit, group = sample_id,
                                    node_groups = node2cluster, normalize = FALSE)
# crumblr wants a plain integer samples-by-cluster matrix, so drop the
# somalign_soft_frequencies class and round the soft counts to integers.
counts <- matrix(as.integer(round(counts)), nrow(counts), ncol(counts),
                 dimnames = dimnames(counts))
cobj <- crumblr::crumblr(counts)   # CLR in cobj$E, precision weights in cobj$weights
```

The most-likely label is unchanged by softening; only the abundance
estimate is smoothed.
