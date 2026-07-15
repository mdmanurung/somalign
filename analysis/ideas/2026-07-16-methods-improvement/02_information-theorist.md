# Information Theorist — 2 ideas

---

# Mutual Information as an Alignment Sharpness Diagnostic and Epsilon Selector

## Persona
**Information Theorist** — the transport plan P is a joint distribution; I(query node; reference node) tells you how much the alignment actually says.

## Motivation

`epsilon` is currently selected by manual inspection of `somalign_sensitivity_grid` — a laborious grid search with no stopping rule. The core problem is that epsilon controls the *entropy* of the transport plan: at large epsilon, P approaches the product measure a⊗b and I(query; reference) → 0 (the plan says nothing about which reference node a query node maps to); at epsilon → 0, P concentrates mass and I → H(query node marginal). The optimal epsilon lives somewhere in between — sharp enough to carry label signal, diffuse enough to regularize against SOM topology noise.

This gives a principled, computable criterion: choose epsilon so that the *mutual information* I(query node; reference node) under the row-normalized plan P̃ = P/sum(P) is maximized subject to a marginal entropy constraint, or equivalently, treat epsilon selection as a rate-distortion problem where distortion is the expected transport cost and rate is I(query; reference). A plot of I vs. expected cost (the "information curve") as epsilon varies gives the analogue of a rate-distortion frontier: the elbow is the operating point where you get most of the alignment information for the least regularization cost. This can replace or complement `somalign_sensitivity_grid` with a single interpretable curve and an automatic elbow-detection rule.

Beyond epsilon selection, I(query; reference) per se is a valuable per-fit scalar diagnostic: it summarizes alignment sharpness in one number (bits) independent of grid size, making fits comparable across datasets, parameter settings, and SOM resolutions. The per-row conditional entropy H(reference node | query node i) flags individual query nodes that remain ambiguously aligned — exactly the OtherT/NK/CD4T boundary nodes where label transfer is unreliable.

## Connection to Existing Code/Data

- `transport$plan` (M×K matrix, stored as `fit$transport_plan`) is the raw joint measure. Row-normalizing it via `.somalign_row_normalize()` gives the conditional P(reference | query). The product marginal a⊗b is `outer(query$node_masses, reference$node_masses) / sum(...)`. Both are immediately available after `.somalign_solve_ot()`.
- `somalign_sensitivity_grid` already loops over epsilon/rho combinations and returns a tidy data frame — the mutual information computation can be added as one extra column per grid point with negligible cost.
- `somalign_diagnostics` is the natural home for per-fit I and per-node H(reference | query = i); these slot alongside `diagnostics$nodes$match_fraction`.
- The `label_transfer$entropy` column already computed in `.somalign_transfer_labels()` is the Shannon entropy of the posterior label distribution per node — a downstream cousin of this idea. Per-node transport entropy is the upstream version, before label aggregation.
- `somalign_sensitivity_grid` returns scalars per grid point; adding I(query; reference) requires only `sum(P̃ * log(P̃ / (a⊗b)))` over nonzero entries of P̃ — three lines of R.

## Approach

1. Add a helper `.somalign_plan_mutual_information(plan, a, b)` in `ot.R` that computes I(query; reference) = Σ_{i,j} P̃_{ij} log(P̃_{ij} / (ã_i b̃_j)) where ã, b̃ are the plan's empirical marginals (not the input masses, to handle unbalanced mass destruction). Also return per-node H_i = -Σ_j P̃_{ij|i} log P̃_{ij|i} as the vector of conditional entropies.

2. Surface these in `somalign_diagnostics`: add `diagnostics$alignment$mutual_information` (scalar, bits) and extend `diagnostics$nodes` with `transport_entropy` (per query node, bits). The per-node entropy column directly identifies ambiguously mapped nodes — flag nodes with H_i > log2(K/4) as "low-confidence transport" in the same way `match_fraction < min_match_fraction` flags low-mass nodes.

3. In `somalign_sensitivity_grid`, compute and return `mutual_information` per grid cell. Add a companion plot function `somalign_plot_information_curve(grid)` that draws the rate-distortion frontier: x-axis = expected transport cost E_P[C], y-axis = I(query; reference), one point per epsilon value, with a simple second-derivative elbow detector (diff of diff on the sorted curve) to annotate the recommended operating point.

4. Expose the elbow as `somalign_select_epsilon(grid)` returning the epsilon value at the elbow, optionally with a `method` argument for "elbow" vs. "entropy_fraction" (choose epsilon where I reaches a user-specified fraction, e.g. 0.9, of its maximum — a softer criterion less sensitive to grid resolution).

## Expected Improvement

- Replaces open-ended grid search with an interpretable stopping rule: users see the information curve, pick the elbow or a target fraction, and have a documented justification.
- Per-node transport entropy in `diagnostics$nodes` gives a finer-grained signal than `match_fraction` alone: a node can have high match fraction but high entropy (mass is transported but spread across many reference nodes — label transfer will be unreliable). This would have directly identified the OtherT/NK/CD4T boundary nodes in the BMV alignment.
- I(query; reference) as a scalar makes alignment quality comparable across datasets and SOM resolutions — essential for benchmarking and for paper-level reporting.

## Feasibility

- **Effort**: Low
- **Fits current architecture**: Yes
- **Methods available**: Standard (Shannon MI, no external dependencies; three lines of R after plan normalization)
- **Key risk**: The plan has many near-zero entries; numerical stability of log(P̃/marginal) requires the same pmax-flooring already used in `.somalign_sinkhorn_kernel`. For the log-domain solver, log P̃ is already computed before `exp()` — a minor code-path optimization available later.

---

# Rate-Distortion Threshold for `outside_reference` via Per-Node Surprisal

## Persona
**Information Theorist** — "outside reference" is an out-of-distribution event; surprisal (negative log-likelihood under the reference SOM's implied density) is a principled OOD score.

## Motivation

The current `outside_reference_distance` flag is triggered when a query cell's distance to its nearest reference SOM node exceeds a node-specific quantile threshold, computed from the reference training data (`reference$distance_quantiles`). This is a reasonable heuristic, but it has two structural weaknesses:

First, it treats all markers equally in Euclidean distance, even though marker-specific batch gain (the BMV CD11c anomaly) or marker-specific noise can dominate the distance for irrelevant reasons. A cell that differs strongly from the reference only on a single marker of known low reliability should not be flagged as outside reference with the same severity as a cell that is uniformly shifted in all markers.

Second, the threshold is node-local but not population-aware: the density of cells around a reference SOM node varies by cell type (tight clusters like naive B cells vs. diffuse clusters like monocytes), but the quantile is computed from raw distances, not from the implied density. A dense cluster has a tight quantile; a diffuse one has a loose quantile. This is partially correct, but the quantile is a scalar that collapses the per-marker structure.

The information-theoretic alternative: model each reference SOM node as a local Gaussian (or diagonal Gaussian, which is free given the reference cell assignments) and compute the **surprisal** — the negative log-likelihood — of each query cell under its assigned reference node's distribution. Surprisal is a proper per-cell scalar in nats (or bits), aggregates over all markers automatically with per-marker variance weighting, and has a chi-squared null distribution (sum of squared z-scores, D markers → chi-squared with D degrees of freedom). This gives a calibrated p-value for "is this cell outside the reference?" with no free threshold parameter: the user chooses a significance level (e.g. 0.01) and the package converts it to a chi-squared critical value.

Per-marker squared z-scores are also decomposable: they tell you *which* markers are responsible for the surprisal, enabling marker-level diagnosis of anomalous cells (e.g. "this cell is outside reference because of CD11c alone").

## Connection to Existing Code/Data

- `reference$distance_quantiles` is computed in `somalign_train_reference` / `_from_som`. The reference cells assigned to each node are already iterated to build these quantiles — the same loop can compute per-node per-marker mean and variance (or equivalently, the diagonal of the within-node scatter matrix) at negligible additional cost during reference construction.
- `reference$codebook` (nodes × features) is the per-node centroid. Per-node variance (nodes × features) is the only additional structure needed.
- `.somalign_project_samples()` already computes the nearest reference node (`unit`) and the Euclidean `distance` for every query cell. Converting from distance to surprisal requires replacing `distance^2` with a per-marker sum of squared z-scores `Σ_d (x_d - μ_{k,d})^2 / σ^2_{k,d}`, which is a one-line vectorized operation once per-node variances are stored.
- `somalign_results()` already outputs `outside_reference_distance` (logical) and `old_som_distance` / `old_som_distance_threshold` (scalars). The surprisal score and p-value slot naturally alongside these — `outside_reference_surprisal` (nats per cell), `outside_reference_pvalue` (from `pchisq`), and `outside_reference_marker_contributions` (per-marker z-scores, as a wide matrix or a top-contributor string).
- The existing `final_status` factor ("inside_reference" / "outside_reference" / "unknown_reference_distance") can be extended with a calibrated variant that uses the p-value threshold instead of the quantile threshold, without breaking the current API.

## Approach

1. During reference construction (in the function that populates `reference$distance_quantiles`), also compute and store `reference$node_var` (nodes × features matrix of per-marker within-node variances, floored at a small epsilon to avoid division by zero for markers with zero within-node variance). Use `reference$node_var <- pmax(node_var, 1e-6)`. This is the only change to the reference object's structure.

2. Add `.somalign_surprisal(x, reference, unit)` in a new or existing utility file: for each cell i assigned to unit k, compute `s_i = Σ_d (x_{i,d} - codebook_{k,d})^2 / node_var_{k,d}` (a chi-squared statistic with D degrees of freedom under a diagonal Gaussian model), the p-value `pchisq(s_i, df = D, lower.tail = FALSE)`, and the per-marker z-scores `z_{i,d} = (x_{i,d} - codebook_{k,d}) / sqrt(node_var_{k,d})`. Vectorized over cells; O(N×D), same as the existing distance computation.

3. Surface in `somalign_results()`: add columns `outside_reference_surprisal` (chi-squared statistic), `outside_reference_pvalue` (calibrated p-value), and `outside_reference_top_marker` (character, name of the marker with the largest |z|-score). Keep the existing `outside_reference_distance` and `final_status` unchanged (backward compatible). Add an optional argument `outside_pvalue_threshold` (default `NULL` to suppress) that, when set, also writes `outside_reference_pvalue_flag` (logical).

4. Add `somalign_plot_surprisal(fit, results)` to `plot.R`: a histogram of per-cell surprisal scores with the chi-squared(D) reference density overlaid, colored by `final_status`. Cells far in the right tail are flagged. A Q-Q plot variant against the chi-squared quantiles immediately shows whether the reference model is well-calibrated (log-linear right tail = Gaussian holds) or heavy-tailed (marker expression is non-Gaussian — expected, and interpretable).

## Expected Improvement

- Replaces an ad hoc quantile threshold with a statistically calibrated OOD score: the p-value threshold corresponds to a known false-positive rate, interpretable to users without cytometry expertise.
- Per-marker surprisal decomposition is a direct diagnostic for marker-specific batch artifacts (like the BMV CD11c anomaly): instead of reporting "N% of cells are outside reference," the package can report "87% of outside-reference cells have their surprisal driven by CD11c alone."
- The diagonal Gaussian model is a strict approximation (SOM nodes are not Gaussian, marker expression is lognormal at best), but the chi-squared reference is still useful as a relative ranking — cells with very high surprisal are genuinely anomalous regardless of model misspecification. The Q-Q plot in Step 4 makes the approximation's quality immediately visible.
- Backward compatible: no existing columns are changed; new columns are additive.

## Feasibility

- **Effort**: Medium
- **Fits current architecture**: Yes, with one small addition to the reference object (`node_var`)
- **Methods available**: Standard (chi-squared distribution from `stats::pchisq`; diagonal Gaussian; no external dependencies)
- **Key risk**: Markers with near-zero within-node variance in the reference (e.g. a binary marker that is uniformly 0 in one cluster) inflate surprisal on query cells that express that marker. The variance floor (`1e-6`) mitigates blow-up, but the p-value will be anti-conservative for such markers. A practical fix: exclude markers whose within-node variance falls below a relative threshold (e.g. < 1% of the global variance) from the surprisal sum and report the effective degrees of freedom.
