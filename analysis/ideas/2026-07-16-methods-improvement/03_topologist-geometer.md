# Topologist / Geometer — 2 ideas

# Laplacian-Regularized Node Shifts: Enforcing a Smooth, Curl-Free Correction Field on the SOM Lattice

## Persona
**Topologist / Geometer** — the SOM grid is a discrete 2D manifold; a physically valid batch correction should be a smooth, low-curl displacement field on it.

## Motivation

The `node_shifts` matrix (M × p, one displacement vector per SOM node) computed in `.somalign_node_shifts()` is currently constructed node-by-node, with zero spatial coupling: each row is the barycentric target minus the query codebook centroid for that node independently. The result is a *displacement field* on the hexagonal SOM lattice — but nothing enforces that adjacent nodes receive similar corrections, or that the field is divergence-/curl-free as expected for a genuine batch shift.

This matters for two reasons rooted in differential geometry:

1. **A true batch effect is a smooth, globally coherent deformation of marker space.** Amplitude gain, reagent lot drift, and laser power fluctuations act on the embedding as a smooth diffeomorphism. The Brenier/Monge theorem, now properly restored by the F2 fix (squared-Euclidean cost in `.somalign_pairwise_distance()`), guarantees that the population-level barycentric map is curl-free in the *continuous* limit — it is the gradient of a convex potential. On the *discrete* SOM grid, this curllessness is only guaranteed in the limit of infinite data and a perfect SOM; in practice, finite-sample OT produces noisy, locally inconsistent shifts. Adjacent nodes that are topologically neighbors (share a hexagonal edge) may receive wildly different correction vectors purely from stochastic mass fluctuation.

2. **Noisy nodes corrupt their cellular neighborhood.** When the corrected codebook is applied in `.somalign_project_pair()`, each cell in `query$sample_unit` inherits the shift of its assigned node. A spuriously large shift on one node artificially displaces thousands of cells in `corrected_matrix`. The current `correction_min_mass` and `match_fraction` gates suppress the worst zero-mass outliers, but do nothing about two neighboring corrected nodes whose shifts differ by, say, 3 z-score units — a geometrically incoherent correction that could push cells across a label boundary.

Regularizing the shift field using the graph Laplacian of the SOM neighbor graph penalizes the squared difference between adjacent shifts, yielding a correction field that is smooth and locally consistent while still fitting the OT-derived targets where data are plentiful.

## Connection to Existing Code/Data

- **`node_shifts`** (`R/fit.R`, `.somalign_node_shifts()`, line 449): the M × p matrix whose rows are the per-node displacements. Currently returned with an attribute `correction_allowed` (logical M-vector). This is exactly the field to be regularized.
- **`shift_transform`** hook (`R/fit.R`, `.somalign_finish_fit()`, lines 143–147): a function slot already present in the pipeline for post-processing `node_shifts`. The subspace projection `function(s) s %*% V %*% t(V)` already uses this hook. A Laplacian smoother can plug in as another `shift_transform`, requiring no architectural change.
- **SOM grid topology**: `kohonen::somgrid` returns a `$pts` matrix of 2D hexagonal coordinates for every node. The neighbor graph (edges between nodes whose `pts` rows are within distance 1 of each other) can be computed in O(M^2) — trivial for typical M = 64–400 nodes. The graph Laplacian L is an M × M sparse matrix derived directly from this neighbor structure.
- **Node masses** (`query$node_masses`): the per-node cell counts, available as weights: high-mass nodes should be fitted more tightly (less regularized) than low-mass nodes.
- **`correction_allowed`**: nodes with `correction_allowed = FALSE` already have `shifts = 0`; in the regularization, these nodes act as zero-valued anchors, pulling neighbors gently back toward zero — a geometrically sensible behavior.

## Approach

1. **Build the SOM neighbor graph and Laplacian.** Extract `query$som$grid$pts` (2D hexagonal coordinates). For a hexagonal grid, neighbors are nodes whose Euclidean coordinate distance equals 1 (up to a small tolerance). Construct the M × M adjacency matrix W (binary), degree matrix D = diag(rowSums(W)), and graph Laplacian L = D − W. This is a one-time O(M^2) computation that costs negligible time for M ≤ 1024.

2. **Formulate the regularized least-squares problem per feature.** For feature dimension d, let s_d ∈ R^M be the raw OT-derived shifts (zeroed for `!correction_allowed` nodes). The regularized shift s*_d solves:

   minimize_{x} ||w ⊙ (x − s_d)||^2 + λ · x^T L x

   where w_i = sqrt(query$node_masses_i) (or w_i = 0 for !correction_allowed). This is a weighted ridge problem with the Laplacian as the regularizer. The closed-form solution is:

   x* = (W_diag + λ L)^{-1} (W_diag s_d)

   where W_diag = diag(w^2). The M × M linear system is solved once per feature dimension p — for M = 256 and p = 40, this is 40 independent small linear solves.

3. **Solve efficiently.** L is sparse symmetric positive semidefinite. The system (W_diag + λ L) is M × M and, for typical M ≤ 1024, can be solved exactly via `solve()` in base R after precomputing the Cholesky factor. For very large SOMs (M > 2000), an iterative conjugate-gradient solver (`stats::conj.grad` shim or sparse matrix via `Matrix` package) is straightforward. Since p features are all solved with the same system matrix, a single Cholesky factorization amortizes across all p right-hand sides.

4. **Wire into the pipeline as a `shift_transform`.** The regularizer is packaged as:

   ```r
   .somalign_laplacian_smooth <- function(shifts, L, lambda, node_masses, correction_allowed) {
     W_diag <- diag(ifelse(correction_allowed, node_masses, 0))
     A <- W_diag + lambda * L
     rhs <- W_diag %*% shifts
     t(solve(A, rhs))  # p × M then transposed to M × p
   }
   ```

   Exposed as a `laplacian_lambda` argument to `somalign_fit()` (default 0, recovering current behavior). Internally, `somalign_fit()` builds L from `query$som$grid` once, then passes the regularizer as `shift_transform` to `.somalign_finish_fit()`.

   The `lambda` parameter has a natural unit: it is in the same scale as the squared marker distances (post-normalization), so `lambda = 0.1`–`1.0` is a natural starting range analogous to `epsilon`.

## Expected Improvement

- **Reduced spurious cross-label corrections.** Adjacent nodes straddling a population boundary (e.g., CD4T/NK border in the BMV dataset) currently receive independent shifts; Laplacian regularization would force these to agree, preventing one node's shift from pulling cells across the label boundary.
- **More stable corrections under query-SOM randomness.** The largest variance source identified in CONTEXT.md is query-SOM training randomness. Laplacian smoothing acts as a geometric prior that absorbs local SOM topology jitter without discarding the global batch signal.
- **Mathematically grounded.** The continuous Monge map for squared-Euclidean cost is the gradient of a convex function, hence curl-free and harmonic except at mass-concentration points. The Laplacian regularizer enforces the discrete analogue of harmonicity on the SOM lattice.
- **Backward compatible.** `laplacian_lambda = 0` (default) reproduces the existing pipeline exactly.

## Feasibility
- **Effort**: Medium
- **Fits current architecture**: Yes — the `shift_transform` hook in `.somalign_finish_fit()` is designed for exactly this; `query$som$grid` exposes the needed geometry
- **Methods available**: Standard — weighted Laplacian smoothing via base-R `solve()`, no new dependencies
- **Key risk**: Choosing `lambda` requires guidance; too large collapses all node shifts to a single global vector (defeating local correction), too small has no effect. A diagnostic plot of smoothed vs. raw shift norms across λ (analogous to `somalign_sensitivity_grid`) would be needed.

---

# Persistent Homology Audit of Node-Shift Topology: Detecting Hole-Inducing Corrections

## Persona
**Topologist / Geometer** — persistent homology tracks how the topology of a point cloud changes across scale; applied to the SOM codebook pre- and post-correction, it reveals whether alignment creates or destroys topological holes (missing populations, disconnected clusters).

## Motivation

A fundamental but currently invisible failure mode of SOM-level batch correction is **topological distortion**: the correction erases a query-only population (a genuine biological cluster becomes a hole) or merges two reference populations (two clusters collapse into one). These events are not caught by any current diagnostic — `outside_reference_distance`, `match_fraction`, and `correction_norm` are all per-node scalars that cannot see multi-node topological events.

Persistent homology (PH) of the codebook point cloud, computed at two scales — pre-correction query codebook and post-correction codebook — provides a scale-free topological fingerprint. Specifically:

- **H_0 (connected components)**: tracks how many clusters exist and when they merge as the neighborhood radius grows. A correction that merges two well-separated clusters will show a 0-cycle in the query that disappears in the corrected codebook.
- **H_1 (loops/holes)**: tracks closed loops in the 1-skeleton of the Vietoris-Rips complex. A correction that maps one population to another while leaving a topological hole may increase 1-cycles.

Comparing the persistence diagrams (birth–death pairs) before and after correction quantifies, for the first time, whether somalign is doing topology-preserving alignment or topology-destroying collapse.

This is geometrically non-trivial: the SOM codebook is a set of M points in R^p (marker space); its topology in marker space may differ from the SOM grid topology (which is always a 2D rectangle or hexagonal mesh). Two SOM nodes that are grid-neighbors may be far apart in marker space if the population landscape is fragmented. It is the marker-space topology that determines biology; the grid topology is an artefact of the SOM training constraint.

## Connection to Existing Code/Data

- **Query and reference codebooks** (nodes × features): the raw point clouds on which PH is computed. Available as `query$codebook` (M × p) and `reference$codebook` (K × p) in reference-scaled space — the same space used by the OT cost `.somalign_pairwise_distance()`.
- **`node_shifts`** (`fit$node_shifts`, M × p): adding these to `query$codebook` yields the corrected codebook. The H_0 and H_1 diagrams can be compared between `query$codebook`, `query$codebook + node_shifts`, and `reference$codebook`.
- **`correction_allowed`** attribute on `node_shifts`: PH should be computed only over the `correction_allowed = TRUE` nodes, since the others are unchanged — zeroing corrected nodes' shifts avoids spurious topology changes from the forced-zero entries.
- **Node masses** (`query$node_masses`): can weight the PH computation — nodes with large mass are more "real" topological features; a filtration weighted by mass creates a more biologically meaningful diagram.
- **`distance_quantiles`** (`reference$distance_quantiles`): the 95th percentile per-node distance threshold provides a natural scale for choosing the PH filtration range. The corrected codebook should have its H_0 persistence diagram compatible with the reference's, meaning populations should not merge at a scale smaller than the within-population spread (captured by `distance_quantiles`).
- **`somalign_diagnostics()`** (`R/diagnostics.R`): the natural home for a new `$topology` diagnostic list containing pre/post PH diagrams, Wasserstein distance between them, and a boolean flag `topology_distorted`.

## Approach

1. **Implement a lightweight Vietoris-Rips PH engine for small point clouds.** For M ≤ 1024 SOM nodes in R^p, the full pairwise distance matrix is already computed by `.somalign_pairwise_distance()` as part of the OT cost. H_0 persistence (equivalent to single-linkage clustering dendrogram) is computed in O(M^2 log M) via sorting the edges and running union-find — entirely in base R with no external C++ dependency. H_1 (loop detection) requires building the Rips complex, which is harder; for the tractable case, restricting to H_0 alone already catches cluster merging/splitting.

   For H_1, a practical approximation: compute the minimum spanning tree (MST) of the pairwise distance graph using Prim's algorithm (O(M^2) in base R), then count the number of back-edges that would form short cycles (length < threshold). This is not true persistent H_1 but catches gross loop creation.

   A stretch option: wrap the `TDA` R package (which uses GUDHI/Dionysus via Rcpp) as an optional dependency checked at runtime — if `TDA` is available, use it for full H_0 and H_1; otherwise fall back to the MST-based H_0 approximation.

2. **Define the pre- and post-correction point clouds.** Three point clouds are compared:
   - `CB_query`: `query$codebook[correction_allowed, ]` (the query nodes that moved)
   - `CB_corrected`: `query$codebook[correction_allowed, ] + node_shifts[correction_allowed, ]`
   - `CB_ref`: `reference$codebook` (the target)

   For each, compute the H_0 persistence diagram (birth = 0, death = merge distance between components). The number of persistence pairs with death > `eps_topology` (a threshold proportional to `max(reference$distance_quantiles["95%"])`) gives the number of "robustly separated populations."

3. **Compute topology change statistics.**
   - `n_components_query`, `n_components_corrected`, `n_components_reference`: count of H_0 pairs with death > threshold.
   - `topology_delta`: `n_components_corrected - n_components_query` (negative = populations merged, positive = populations split by correction).
   - `bottleneck_distance_H0`: Wasserstein-infinity (bottleneck) distance between the H_0 diagrams of `CB_corrected` and `CB_ref`. This measures how well the corrected codebook topology matches the reference topology.
   - A flag `topology_warning` if `|topology_delta| > 0` or `bottleneck_distance_H0 > user_threshold`.

4. **Expose as a `$topology` slot in `somalign_diagnostics()`.** The function is already called at the end of `.somalign_finish_fit()`. Add:

   ```r
   topology = .somalign_topology_audit(
     query$codebook, node_shifts, reference$codebook,
     attr(node_shifts, "correction_allowed"),
     reference$distance_quantiles
   )
   ```

   The returned list contains the persistence diagrams (as data frames of birth/death pairs), the three component counts, the bottleneck distance, and the boolean warning. A companion `somalign_plot_topology(fit)` function renders the three persistence diagrams side-by-side.

## Expected Improvement

- **First-ever topology diagnostic for SOM-level alignment**: currently somalign has no way to detect whether a correction merged two populations or erased one. This fills that gap directly.
- **Early warning for over-correction.** If `laplacian_lambda = 0` and `epsilon` is large (diffuse transport), the corrected codebook can have fewer distinct populations than the reference. The topology audit makes this visible before downstream label-transfer analysis.
- **Interpretable output.** The number of "robustly separated populations" (H_0 components above threshold) is directly biological: CD4T, NK, and OtherT are separate components in the reference; if the corrected codebook collapses NK and OtherT into one, `n_components_corrected` decreases. This would have been a useful early diagnostic for the hard CD4T/NK boundary problem described in CONTEXT.md.
- **Integration with `somalign_sensitivity_grid()`**: recording `bottleneck_distance_H0` across the epsilon/rho grid reveals which parameter regimes are topology-safe.

## Feasibility
- **Effort**: Medium (H_0 only via MST/union-find, no external deps) or High (full H_0 + H_1 with `TDA` wrapper)
- **Fits current architecture**: Yes — plugs into the existing `somalign_diagnostics()` output list; no changes to the fit pipeline itself
- **Methods available**: Standard for H_0 (union-find/MST, base R); Research-grade for H_1 without external dependencies; Standard with `TDA` package as optional dep
- **Key risk**: For p > 20 markers, the pairwise M × M distance is already computed for the OT cost but PH in high-dimensional spaces can be dominated by the curse of dimensionality — topological holes that are geometrically real in 2D PCA may be invisible in 40D marker space. Users may need to run PH on a low-dimensional projection (first 5 PCs of the codebook) rather than the full marker space, adding a parameter choice.
