# Representation Learning Specialist — 2 ideas

---

# Learned Mahalanobis Metric for Batch-Aware OT Cost

## Persona
**Representation Learning Specialist** — fix Euclidean blindness: learn a marker-space metric where batch directions are cheap and biological directions are expensive to transport.

## Motivation

The current OT cost in `.somalign_align_transport` (R/fit.R, line 210) is plain squared Euclidean distance between query and reference codebook rows, normalised by a scalar median (`.somalign_prepare_cost`, line 185). This treats every marker dimension identically: a one-unit shift in CD3 costs the same as a one-unit shift in CD11c, regardless of whether that direction is dominated by batch or biology.

This is precisely backwards for batch correction. Batch effects in CyTOF are marker-specific gain/offset drifts. The anchor displacement matrix D (computed in `.somalign_batch_subspace`, R/anchored.R, line 248) already tells us which directions in marker space are "batch" (high variance in D) versus "biology" (low variance in D, orthogonal to batch subspace V). A Mahalanobis metric Σ^{-1} built so that batch directions have low weight (cheap to transport) and biological directions have high weight (expensive to transport) would concentrate OT mass on corrections that are credibly batch-driven, while penalising transport that would require crossing genuine biological differences.

The real-data context reinforces this: the BMV vs pilot alignment showed a gain/amplitude effect on lineage markers. Some markers carry the batch signal strongly; others are nearly clean. A diagonal Mahalanobis (per-marker weights) would naturally up-weight the clean markers as anchors for transport geometry and down-weight the noisy ones, reducing over-correction on biology-carrying markers.

## Connection to Existing Code/Data

- **Cost computation**: `.somalign_pairwise_distance` (called at fit.R line 210) returns M × K squared Euclidean distances. The fix is a single drop-in replacement that computes the Mahalanobis cost `(q_i - r_j)^T Σ^{-1} (q_i - r_j)` where Σ^{-1} is a learned p × p (or diagonal p-vector) precision matrix.
- **Anchor displacement matrix D**: `.somalign_batch_subspace` (anchored.R, line 248) computes `anchor_old_scaled - anchor_new_scaled`. The batch subspace V (a p × rank matrix) is already returned as `fit$anchors$batch_subspace$V`. The full SVD of D is already available as input to `.somalign_subspace_svd`.
- **Reference center/scale**: `reference$center` and `reference$scale` define the standardised space where the cost is already computed; the Mahalanobis precision operates in this same space, so no coordinate change is needed.
- **No-anchor fallback**: for `somalign_fit` (no anchors), a diagonal Σ^{-1} estimated from `reference$codebook` column variances (or a user-supplied weight vector) suffices, analogous to `diagonal_boost` which already tweaks cost per-node (fit.R line 191).

## Approach

1. **Estimate precision from anchors.** Given D (n_anchors × p), compute a diagonal precision Σ^{-1} = diag(w) where `w_j = 1 / (var(D[,j]) + delta)` with delta a ridge term. Dimensions where anchor displacement variance is large (batch-driven) get small weight (cheap transport); dimensions where D is near-zero (biological) get large weight (expensive transport). Optionally use the full SVD to build a low-rank + identity precision: `Σ^{-1} = I + (lambda - 1) V V^T` where lambda < 1 down-weights the batch subspace and lambda > 1 up-weights the biological complement.
2. **Implement `.somalign_mahalanobis_cost(query_codebook, reference_codebook, precision_diag)`** that computes `rowSums(sweep((q_i - r_j)^2, 2, w, "*"))` for diagonal case — no matrix inversion needed, just element-wise scaling before the sum-of-squares. For the full M × K matrix this is equivalent to computing distances in a whitened space: `scale each column j of the codebook by sqrt(w_j)`, then call the existing `.somalign_pairwise_distance`. This is a one-liner on top of existing infrastructure.
3. **Add `metric` argument to `somalign_fit_anchored`.** Options: `"euclidean"` (default, current behaviour), `"mahalanobis_diagonal"` (auto-estimated from anchor D), `"mahalanobis_subspace"` (low-rank + identity from batch subspace V), `"user"` (user-supplied weight vector). The precision is estimated once from anchors before the OT solve and stored in `fit$anchors$precision_diag` for inspection.
4. **Expose `marker_weights` parameter on `somalign_fit` (no anchors)** for the anchor-free case. Default NULL (Euclidean). When supplied as a named numeric vector (one per feature), apply diagonal Mahalanobis. This lets users who know their noisy markers supply informed weights without anchors.

## Expected Improvement

- Transport mass concentrates on biologically concordant node pairs rather than being diffused by noisy/batch-driven markers.
- Barycentric correction shifts become smaller in biological directions (the expensive dimensions), reducing over-correction artefacts seen in the BMV alignment's OtherT/NK/CD4T boundary region.
- The learned metric is interpretable: `fit$anchors$precision_diag` directly shows which markers are batch-driver vs biology-carrier, a useful diagnostic in its own right.
- For the hard case where anchors are unavailable, `marker_weights` provides a principled lever currently absent (only `diagonal_boost`, which modulates the nearest-reference-node pair globally, not per marker dimension).

## Feasibility

- **Effort**: Medium
- **Fits current architecture**: Yes — `.somalign_pairwise_distance` is a one-function change; the precision matrix is estimated before the OT solve and passed through `.somalign_prepare_cost`; all downstream code (Sinkhorn, barycentric map, diagnostics) is unchanged.
- **Methods available**: Standard — diagonal Mahalanobis is textbook; the low-rank variant uses the V already computed by `.somalign_batch_subspace`.
- **Key risk**: With small n_anchors (< p), variance estimates of individual markers in D are noisy. The ridge term delta is a new hyperparameter. A cross-validation scheme (hold out anchor pairs) could calibrate it, but for a first implementation a fixed delta = 1e-2 in reference-scaled space (where codebook entries are O(1)) is defensible. Also: the diagonal variant changes the OT objective, which may shift the transport plan in ways users find surprising if they have tuned epsilon empirically — clear documentation of the interaction is needed.

---

# Anchor-Free Batch Subspace via Contrastive Multi-Batch Invariance

## Persona
**Representation Learning Specialist** — disentangle batch from biology without anchors by exploiting invariance: what is shared across many query batches is biology; what differs is batch.

## Motivation

The batch subspace V in `somalign_fit_anchored` (anchored.R line 222) requires the same biological sample to be measured in both batches — an anchor. In many cytometry studies no such paired measurements exist. The two-pass diagnostic in `somalign_fit_two_pass` (fit.R line 647) gives a read-only batch subspace estimated from OT-derived node shifts, but explicitly warns it "may conflate batch effects with biology" and is not used for correction.

A principled anchor-free approach: if a user has K ≥ 2 query batches all being aligned to the same reference, then the directions in codebook space that differ most across batches (high between-batch variance of node shifts) are likely batch directions, while directions with low between-batch variance but large departure from zero are biology. This is a contrastive decomposition: shared signal (consistent shift across batches) = biology change from reference; variable signal (batch-specific shift) = batch nuisance.

This extends somalign from a 1-vs-1 alignment paradigm to a K-vs-1 multi-batch paradigm where batches collectively supervise each other's batch-subspace estimation — without any anchor samples.

## Connection to Existing Code/Data

- **Node shifts per batch**: `.somalign_node_shifts` (fit.R line 449) returns an M × p matrix of per-node corrections (barycentric target minus query codebook). For K batches aligned separately to the same reference, we get K such matrices `NS_1, ..., NS_K`, each M × p.
- **Two-pass batch subspace**: `.somalign_subspace_svd` is already implemented and called at fit.R line 648 on node shifts from pass 1. The function takes a weighted matrix and a variance threshold — exactly what is needed.
- **`somalign_fit_two_pass`**: the pass-1 node shifts `ns1` (fit.R line 593) are already the best available anchor-free proxy for the batch direction. The multi-batch version generalises this from one pass to K batches.
- **`fit$node_shifts`**: stored in every `somalign_fit` object (fit.R line 357), making it easy to collect these matrices from K independent fits before the correction step.

## Approach

1. **New function `somalign_fit_multibatch(query_list, reference, ...)`** that accepts a named list of K `somalign_query` objects (one per batch), runs a first-pass OT fit for each against the shared reference (using large `epsilon_global`), and collects the K pass-1 node-shift matrices `NS_1, ..., NS_K` (each M × p, restricted to nodes that are correction-allowed in at least one batch).
2. **Estimate the batch subspace** from the between-batch variation of node shifts. Stack the K shift matrices and decompose: compute the node-mass-weighted column covariance of `NS_k - NS_mean` across batches, where `NS_mean` is the mass-weighted mean shift (the "shared" signal). SVD of this between-batch difference matrix gives V_batch (directions that vary across batches = batch nuisance). The complement of V_batch in the node-shift space contains the consistent shift = biology-from-reference.
3. **Apply subspace correction to each batch independently**: run a second-pass OT for each query at smaller `epsilon_local`, then project node shifts onto V_batch (same mechanics as `shift_fn <- function(s) s %*% V %*% t(V)` in anchored.R line 224). This removes the batch-variable component and retains only the biologically consistent correction.
4. **Return a `somalign_multibatch_fit` object**: a named list of K `somalign_fit` objects sharing the same `batch_subspace$V`, plus a summary slot with between-batch variance explained per marker and rank of V_batch. Each per-batch fit is fully compatible with `somalign_results()`.

## Expected Improvement

- Enables principled subspace-constrained correction without anchor samples — the current blocker for using `correction = "subspace"` in the common case where paired measurements are unavailable.
- The shared-biology vs batch-variable decomposition is interpretable and can serve as a quality metric: if the between-batch variance is near zero (small V_batch rank), the batches are already well-aligned and little correction is warranted; if it is large and low-rank, a clean batch direction exists.
- The multi-batch estimate of V_batch is more robust than the two-pass single-batch diagnostic: averaging K between-batch differences suppresses node-shift noise that would otherwise masquerade as batch directions.
- Directly addresses the known limitation "full-space barycentric correction erases query-only biology" for the anchor-free case.

## Feasibility

- **Effort**: High
- **Fits current architecture**: Needs refactor — a new top-level function is required; the internals (`.somalign_node_shifts`, `.somalign_subspace_svd`, `shift_fn` pattern from anchored.R) are fully reusable with no changes. The main engineering effort is the multi-batch orchestration loop and the between-batch covariance estimation step.
- **Methods available**: Standard — weighted PCA / SVD of between-batch shift differences is textbook. The `.somalign_subspace_svd` helper already encapsulates this pattern.
- **Key risk**: The method requires K ≥ 2 batches all measured on the same panel against the same reference — a reasonable assumption in longitudinal CyTOF studies (the BMV use-case has multiple time-point batches) but not universally available. With K = 2 the between-batch covariance is rank-1 and V_batch is a single vector, which may be too coarse. A minimum K = 3–5 is advisable; a warning should be emitted for K < 3. Additionally, if one batch has genuine biological novelty (new cell populations absent in others), that signal may project onto V_batch and be erroneously suppressed; a node-level flag of "novel nodes" (high `outside_reference_distance`) should gate those nodes out of the V_batch estimation.
