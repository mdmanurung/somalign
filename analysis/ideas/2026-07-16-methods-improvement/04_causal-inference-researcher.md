# Causal Inference Researcher — 2 ideas

---

# Sensitivity Analysis for Unmeasured Batch Confounding (Rosenbaum-Style Bounds on the Correction)

## Persona
**Causal Inference Researcher** — the barycentric correction estimates a counterfactual ("what would this cell look like in the reference batch?"), so we should quantify how sensitive that estimate is to violations of the identifying assumptions.

## Motivation

The identifiability argument for somalign's batch correction rests on a causal graph: **batch → expression ← biology**, where batch is a confounder between the query and reference distributions. The anchors (same biological sample in both batches) function as a matched-pair negative control that identifies the batch direction D = anchor_old − anchor_new without biological confounding — but only under the exclusion restriction that anchor measurements carry *only* batch variation, not biology-×-batch interactions or anchor-specific sample artifacts.

If this restriction is violated — e.g., anchor samples are not representative of the full biology, or there is an unmeasured marker-specific gain that differs between populations — the subspace V estimated from SVD(D) is a biased instrument, and every node shift projected onto V inherits that bias. Currently the package gives no way to quantify how large this violation can be before the correction is meaningless. Sensitivity analysis (in the Rosenbaum sense: "how strong must the unmeasured confounding be to explain away the effect?") turns this from a binary "we trust the anchors" assumption into a calibrated, reportable quantity.

Concretely: the package already exposes `fit$anchors$batch_subspace` (V, rank, variance_explained) and the unconstrained node shifts `fit$node_shifts`. A practitioner can ask: if the true batch subspace deviates from V by angle θ (in the Grassmannian sense), how much does the corrected codebook coordinate of each node change? For small θ the correction is stable; for large θ it is meaningless. Reporting the "tipping-angle" — the smallest angular perturbation to V that would reverse the sign of the correction for the median node — gives a decision-relevant sensitivity summary.

## Connection to Existing Code/Data

- `fit$anchors$batch_subspace$V` (p × rank): the estimated batch subspace from `.somalign_batch_subspace()` → `.somalign_subspace_svd()`. V is the instrument.
- `fit$node_shifts` (M × p): raw barycentric shifts from `.somalign_node_shifts()` before projection onto V. The projected shift for node i is `shift_transform(node_shifts[i,]) = node_shifts[i,] %*% V %*% t(V)`.
- `fit$anchors$n_anchors`, `fit$query$node_masses`: degrees of freedom and weighting for bootstrap resampling of D.
- `somalign_sensitivity_grid()` already exists as a precedent for a "grid of perturbations → diagnostic summary" pattern; this idea extends it to the subspace direction rather than epsilon/rho.

## Approach

1. **Define the sensitivity parameter.** For a perturbation matrix ΔV (p × rank) with small Frobenius norm δ, form a perturbed subspace V(δ) = orth(V + ΔV) (orthonormalised). The perturbed correction for node i is `node_shifts[i,] %*% V(δ) %*% t(V(δ))`. The sensitivity of the correction is dC/dδ evaluated in the worst-case direction ΔV.

2. **Bootstrap the anchor displacement matrix D.** Resample rows of D = anchor_old_scaled − anchor_new_scaled (already computed inside `.somalign_batch_subspace()`) with replacement B times (default B = 200). Recompute SVD(D_b) and the resulting V_b for each bootstrap replicate. This gives an empirical distribution over V that reflects anchor sampling variance without requiring analytical derivatives.

3. **Compute bootstrap node-shift distributions.** For each replicate b, apply `node_shifts %*% V_b %*% t(V_b)` and compare to the point estimate. Store per-node 2.5th/97.5th percentiles of the correction norm and of each feature's corrected coordinate. The width of these intervals is the empirical sensitivity to anchor sampling.

4. **Report tipping angles and a summary diagnostic.** Compute the principal angle between V and each V_b (using `svd(t(V) %*% V_b)$d` → `acos(clamp(d, -1, 1))`). Plot correction-norm CI width vs. principal angle. Add a `$sensitivity` slot to the `somalign_anchored_fit` object containing: `bootstrap_V` (array p × rank × B), `node_correction_ci` (M × 2), `median_tipping_angle_deg`, and `anchor_leverage` (the Cook's-distance analog: which anchor pairs most influence V). Expose a new S3 method `somalign_subspace_sensitivity(fit, B = 200, seed = 1L)`.

## Expected Improvement

Users gain a principled, reportable answer to "how much do we trust the batch correction?" anchored in the same causal logic as the subspace mode. High sensitivity (wide CIs, small tipping angle) flags cases where more or better-chosen anchors are needed before publication. Low sensitivity provides positive evidence that the correction is robust to the exclusion restriction. This is directly analogous to the E-value in epidemiological sensitivity analysis, adapted to the subspace-estimation context.

## Feasibility

- **Effort**: Medium
- **Fits current architecture**: Yes — adds one new exported function `somalign_subspace_sensitivity()` and a new `$sensitivity` slot; does not modify any existing path; follows the pattern of `somalign_sensitivity_grid()`.
- **Methods available**: Standard — bootstrap resampling and principal-angle computation are textbook linear algebra; no new dependencies needed beyond base R `svd()`.
- **Key risk**: With small anchor counts (n_anchors < rank × 5) the bootstrap distribution degenerates and the sensitivity interval is uninformative; the function should warn and recommend a minimum anchor-to-rank ratio.

---

# Anchor Exclusion-Restriction Test via Conditional Independence of Residuals

## Persona
**Causal Inference Researcher** — anchors are a valid instrument only if the residual (biology-direction) component of the anchor displacement is independent of the batch label; this is the exclusion restriction, and it is testable.

## Motivation

In instrumental-variables (IV) language, the anchors identify the batch effect under two conditions: (1) relevance — the anchor displacement D = anchor_old − anchor_new has substantial variance in the batch direction (guaranteed when variance_explained is high), and (2) exclusion restriction — the component of D orthogonal to the batch subspace V is zero in expectation (i.e., anchors carry no systematic biology-direction displacement). Condition (2) is the fragile one: if the anchor samples were, for example, preferentially collected at a particular stage of sample preparation, they may carry a biology-×-batch interaction that inflates D in a biology direction.

The exclusion restriction is currently assumed but untested. Yet it is testable: under the null (valid instrument), the projection of D onto the orthogonal complement V_perp = I − V %*% t(V) should be a zero-mean matrix with no structure — specifically, it should not correlate with any external covariate of the anchor samples (time of measurement, operator, sample concentration), and its singular values should not exceed what is expected under isotropic noise of the same variance. A simple test: compare the leading singular value of `D %*% V_perp` to a permutation null (shuffle anchor identities, recompute D_perm %*% V_perp, take leading SV). A p-value below 0.05 means the residual has detectable structure that the batch subspace did not absorb — a red flag for exclusion-restriction violation.

This is directly analogous to the overidentification test (Sargan–Hansen J-test) in classical IV regression, adapted to the high-dimensional subspace setting.

## Connection to Existing Code/Data

- `anchors_scaled$anchor_old_scaled` and `anchors_scaled$anchor_new_scaled` (n_anchors × p), available inside `.somalign_anchored_dispatch()` after `.somalign_validate_anchors()`, passed to `.somalign_batch_subspace()`. Currently these are discarded after V is computed; they need to be stored or re-exposed.
- `fit$anchors$batch_subspace$V` (p × rank): the instrument. The orthogonal complement projector is `diag(p) - V %*% t(V)`, computable from V alone with no new data.
- `fit$anchors$n_anchors`: degrees of freedom for the permutation test.
- The `somalign_diagnostics()` function provides precedent for attaching structured test results to a fit object.

## Approach

1. **Store anchor displacements in the fit object.** Extend `fit$anchors` to include `D_scaled` (n_anchors × p), the scaled displacement matrix. This is a small (n_anchors × p) matrix — for typical CyTOF data n_anchors < 200, p < 50 — so storage overhead is negligible. Modify `.somalign_anchored_dispatch()` to pass `D_scaled = anchors_scaled$anchor_old_scaled - anchors_scaled$anchor_new_scaled` through to the returned fit.

2. **Compute the residual matrix.** In a new exported function `somalign_exclusion_test(fit, n_perm = 999L, seed = 1L)`, form `R = D_scaled %*% (diag(p) - V %*% t(V))` (n_anchors × p), then `sv_obs = svd(R, nu = 0, nv = 0)$d[1]` (leading singular value of the residual).

3. **Permutation null.** Shuffle the rows of `anchor_old_scaled` (breaking the within-anchor-pair correspondence while preserving marginal distributions) n_perm times; for each permutation recompute `D_perm = anchor_old_perm - anchor_new_scaled`, `R_perm = D_perm %*% V_perp`, and record `sv_perm[b] = svd(R_perm)$d[1]`. The permutation p-value is `mean(sv_perm >= sv_obs)`.

4. **Return a structured result and integrate into diagnostics.** Return a list: `sv_observed`, `sv_null_quantiles` (2.5/50/97.5th percentiles of sv_perm), `p_value`, `rank_used` (the subspace rank), `verdict` ("pass"/"warn"/"fail" thresholded at p > 0.1 / 0.05 / 0.01). Optionally: a per-feature decomposition of the residual norm to identify which markers drive the violation. Attach the result to `fit$anchors$exclusion_test` when the user calls `somalign_exclusion_test()`, and print a one-line summary in the `somalign_fit_anchored()` output when `correction = "subspace"` or `"both"` and n_anchors is large enough (≥ 10) to make the test informative.

## Expected Improvement

Users get the first formal validity check for the anchor instrument — the causal analogue of a model diagnostic rather than just a goodness-of-fit metric. A failing test (p < 0.05) tells the user that the subspace mode is removing biology as well as batch, and they should either increase anchor diversity, switch to `correction = "cost_bonus"`, or investigate which markers are driving the residual. A passing test provides positive evidence that the IV assumptions hold, strengthening the scientific credibility of the batch correction for publication. This fills the most important gap in the current causal argument for the subspace correction.

## Feasibility

- **Effort**: Low
- **Fits current architecture**: Yes — one new exported function, one small addition to the `$anchors` slot (storing D_scaled), no changes to any existing computational path.
- **Methods available**: Standard — permutation test with SVD is base-R and runs in seconds even for B = 999 with small n_anchors.
- **Key risk**: With very few anchors (n_anchors < rank + 2) the test has essentially no power; the function should emit a warning and suggest a minimum of 3 × rank anchor pairs for the permutation null to be meaningful.
