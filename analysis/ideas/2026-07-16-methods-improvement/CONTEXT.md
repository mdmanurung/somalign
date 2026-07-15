# somalign — package context for ideation

**Goal of this session:** generate ideas to *improve the somalign R package* (new
features, algorithms, diagnostics, robustness, theory) — NOT to analyze a dataset.
Every idea must be grounded in the actual architecture and data structures below.

## What somalign is
An R package (Bioconductor-style, v0.99.1) that aligns a **query** cytometry/CyTOF
batch onto a **fixed reference** batch at the level of **SOM codebooks**, using
**codebook-level unbalanced entropic optimal transport (UOT)**. The direct
reference projection is the conservative primary result; transport-corrected
projections are auxiliary (annotation/visualization). Pure-R Sinkhorn; depends on
`kohonen`, `ggplot2`, `stats`; optional `BiocParallel`.

## Pipeline
1. **Reference SOM** (`somalign_reference` / `_train_reference` / `_from_som`):
   train a Kohonen SOM on reference data (winsorized + standardized via
   `center`/`scale`). Store: codebook (nodes × features) in *reference-scaled*
   space, node masses, per-node distance quantiles (thresholds for "outside
   reference"), optional per-node label probabilities (from a supervised xyf SOM
   second layer).
2. **Query** (`somalign_query` / `_from_som`): assign query cells to a query SOM;
   scale the query codebook into reference-scaled space using `reference$center/scale`.
   Node masses = assignment counts.
3. **OT** (`.somalign_align_transport`, `R/ot.R`): cost = **squared** Euclidean
   distance between query and reference codebooks, normalized by median positive
   entry. Solve UOT (Sinkhorn) with entropic reg `epsilon` and KL marginal
   relaxations `rho_query`, `rho_ref`. Two solvers: standard scaling
   (`K = exp(-C/eps)`) and numerically stable `log_domain`.
4. **Node shifts** (`.somalign_node_shifts`, `R/fit.R`): barycentric map — each
   query node shifted toward the mass-weighted mean of reference nodes under the
   row-normalized transport plan. `correction_allowed` gates which nodes move.
5. **Correction + outputs** (`somalign_results`): shift query codebook, reproject
   cells. Columns: `old_som_unit` (direct), `corrected_som_unit`,
   `outside_reference_distance`, `correction_norm`, `match_mass_ratio`,
   `transferred_label` (+ `confidence`, `accepted`, and new `_second` /
   `_second_confidence` / `_margin` for triaging close argmax ties).

## Key features
- **Anchor-guided** (`somalign_fit_anchored`): anchors = the *same* biological
  sample measured in both batches. `correction` modes: `"cost_bonus"` (anchor
  counts bias OT cost toward anchor-supported node pairs), `"subspace"`
  (CellANOVA-inspired: SVD of anchor displacements → batch subspace V; restrict
  node shifts to V so biology orthogonal to V is preserved), `"both"`.
- **Two-pass** (`somalign_fit_two_pass`): global shift (large eps) + residual
  (small eps); read-only batch-subspace diagnostic.
- **Diagnostics**: `somalign_diagnostics`, `somalign_sensitivity_grid`
  (epsilon/rho grid), `somalign_som_stability` (across query-SOM seeds),
  `somalign_check_codebook_alignment` (range overlap, centroid drift, cost
  coverage before per-cell work).
- **Plots** (`R/plot.R`): match fraction, label confusion, etc.

## Data structures you can build on
Reference/query codebooks (nodes × features), node masses (a, b), per-node
distance quantiles, the **M×K transport plan** P, node shifts (M × features),
label-probability matrices (nodes × classes), per-cell results, anchor
displacement matrices D (n_anchors × features), batch subspace V (features × rank).

## Recent work (this session)
- Fixed 7 scientific-correctness issues, notably: **ANCHOR-001** (anchor codebook
  projections were swapped — critical), **F2** (OT cost switched to squared
  Euclidean, restoring the Brenier optimal-transport-map property of the
  barycentric correction), F3 (all-underflow Sinkhorn-row warning), Inf-threshold
  handling. All tests pass; pushed.
- Real use: aligning a 39.8M-cell "BMV" CyTOF batch onto a "pilot" reference.
  Found a gain/amplitude batch effect on lineage markers, fixed by
  quantile-normalization (divide by the query's own 99.9th percentile before
  applying the reference's center/scale). Hard residual: OtherT/NK/CD4T label
  boundary; an anomalous CD11c signal in query CD4T cells.

## Known limitations / open problems (fertile ground)
- Full-space barycentric correction **erases query-only biology** (motivated the
  subspace mode, but that needs anchors).
- Label-transfer argmax is **brittle at close ties**; `_second`/`_margin` are
  node-level, only a partial fix.
- **No uncertainty quantification** on node shifts or label transfer.
- **Query-SOM training randomness** is the largest uncontrolled variance source
  (`som_stability` measures it but doesn't fix it).
- `epsilon`/`rho` selection is **manual** (sensitivity grid, no principled rule).
- Cost is plain Euclidean in marker space — **no learned/Mahalanobis metric**, no
  per-marker weighting, no handling of marker-specific batch gain.
- **Single global transport**; no local, hierarchical, or per-population structure.
- Underflow at small `epsilon` (log-domain solver mitigates but is slower;
  pure-R, no C++).
- Correction is applied to codebooks then cells re-projected — **no per-cell
  transport**, so within-node structure is lost.
