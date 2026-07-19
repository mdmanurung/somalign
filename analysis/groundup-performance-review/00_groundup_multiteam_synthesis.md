# Ground-Up Multi-Team Synthesis

Date: 2026-07-19

## Scope

This document records a multi-team review of SOMalign package performance from
first principles. It integrates four independent read-only reviews:

1. Mathematical/statistical formulation.
2. Biological label-transfer behavior, focused on the OtherT/T/NK boundary.
3. Validation and publication readiness.
4. Runtime and API engineering.

The goal is to identify the next package improvements that should raise label
transfer precision/recall, strengthen repeat-sample frequency reproducibility,
and make release claims easier to defend.

## First-Principles Model

SOMalign is currently a node-level alignment system:

1. Reference and query cells are scaled into the reference feature space.
2. Each dataset is summarized by SOM nodes: codebook vectors and node masses.
3. Query nodes are aligned to reference nodes with unbalanced entropic OT.
4. The transport plan is row-normalized into a query-to-reference
   correspondence matrix.
5. Query-node label probabilities are computed as:

```text
Pr(label | query node i) = sum_j correspondence[i, j] * Pr(label | reference node j)
```

6. Cells inherit labels from their assigned query SOM node.
7. Soft abundance is separate: each query cell is smoothed over nearby reference
   nodes and aggregated to sample-level frequencies.

This gives the package a clean and scalable contract, but it also explains the
main failure modes. If a query node mixes several T/NK states, its transferred
label posterior is diluted before any cell-level decision is made. If `epsilon`
is too low, recall can improve while precision drops. If `epsilon` is too high,
precision among accepted cells can look good while true class recall and
coverage fall.

## Current Evidence

### Repeat-Sample Frequencies

The strongest current empirical result is repeat-sample abundance
reproducibility from the SOMalign example notebook metadata:

- Repeat metadata: 20 matched pilot-vs-BMV sample pairs.
- Hard nearest-node metacluster CLR weighted-r median: 0.9322.
- Soft kNN metacluster CLR weighted-r median: 0.9566.
- Soft kNN improves all 20/20 matched repeats.
- Lineage-level correlations are near ceiling: hard 0.9873, soft 0.9938.
- Split-half sampling ceiling is approximately 0.9985.

Interpretation: soft abundance should be treated as the preferred estimator for
cluster frequency reporting. This is a frequency result, not proof that hard
cell-level label transfer improved.

### OtherT Label Transfer

Current Batch1-to-Batch2 labelled BMV pilot experiment at 30k cells, 10x10 SOM,
`rlen = 20`:

| Setting | OtherT precision | OtherT recall | OtherT F1 |
|---|---:|---:|---:|
| Baseline `epsilon = 0.1` | 0.9653 | 0.7236 | 0.8272 |
| Global `epsilon = 0.02` | 0.9140 | 0.8632 | 0.8879 |
| T-marker weights, `epsilon = 0.05` | 0.9306 | 0.8528 | 0.8900 |
| `epsilon = 0.05`, `rho = 2.2` | 0.9542 | 0.8176 | 0.8807 |
| High `epsilon = 0.2` | 0.9994 | 0.6568 | lower recall |

Interpretation: the package can move along the OtherT precision/recall curve,
but no deployable single global setting has robustly improved both precision and
recall over the baseline. The most promising signal is marker-aware T/NK
geometry plus targeted rescue, not another global epsilon default.

### Publication Readiness

Current package check status is clean in `somalign.Rcheck/00check.log`. The
validation surface is already broad, with tests and helpers for cross-validation,
label metrics, calibration, tuning, soft labels, sensitivity, and anchors.

The publication-readiness gap is that the strongest quantitative claims still
depend on ad hoc analysis artifacts. They need a machine-readable validation
matrix with raw predictions, confidence intervals, source artifact hashes, and
explicit pass/fail gates.

### Runtime

The runtime bottleneck is not OT. The package-owned hot path is dense
nearest-code projection: an `n_cells x n_nodes` distance matrix is materialized.
SOM training also dominates end-to-end runtime, but that is primarily delegated
to `kohonen::som()`.

High-priority engineering improvements are chunked constructors, reuse of
`som$unit.classif`, compact object storage, streaming corrected projection, and
explicit `BPPARAM`/control objects for memory-aware sweeps.

## Consensus Diagnosis

The package is scientifically coherent and close to release-quality, but its
performance story needs three upgrades:

1. Make uncertainty and calibration first-class, rather than relying on accepted
   hard labels.
2. Treat difficult residual labels such as OtherT as local biological boundary
   problems, not as a single global OT tuning problem.
3. Turn the current analyses into automated release gates with provenance.

## Ranked Improvement Program

### 1. Calibrated Multi-Objective Tuning

Tune `epsilon`, `rho_query`, `rho_ref`, confidence thresholds, and margin
thresholds against full metrics: `accuracy_all`, full macro-F1, coverage,
ECE/Brier, and per-class precision/recall. Accepted-cell MCC must never be the
primary objective alone.

Success criterion: selected settings improve or preserve full macro-F1,
coverage, and calibration across held-out folds and Batch1-to-Batch2 transfer.

### 2. Strict Biological T/NK Router

Replace the broad top-or-second-label T/NK router with a high-specificity router.
Use deployable evidence only: CD3/CD7 or NK-marker support, negative B/myeloid
markers, global margin, entropy, outside-reference diagnostics, direct reference
projection, and soft probabilities.

Success criterion: OtherT recall improves while B-cell and Mono/DC recall do not
collapse, and non-T/NK false routing is sharply reduced from the current broad
router.

### 3. Targeted OtherT Rescue

Start from the high-precision baseline and selectively rescue OtherT candidates
from low-epsilon and T-marker-weighted projections. Rescue only candidates with
strong OtherT posterior, T-cell evidence, ambiguous global margin or abstention,
and no clean CD4/CD8/NK signature.

Success criterion: OtherT recall increases at fixed precision thresholds
(`>= 0.95` and `>= 0.97`) without lowering full macro-F1.

### 4. Shared Soft Posterior for Labels and Frequencies

Unify label transfer and soft abundance through a per-cell posterior:

```text
cell -> query-node softness -> OT correspondence -> reference label probabilities
```

Optionally blend this with direct reference-kNN posteriors. The output should
support hard labels, calibrated abstention, and frequency estimates from the
same probability object.

Success criterion: repeat-sample CLR correlations meet or exceed current soft
kNN performance, while hard-label macro-F1 and OtherT F1 do not regress.

### 5. Reference-Only Feature/Metric Learning

Learn diagonal feature weights from reference labels, class-boundary separation,
anchor consistency, or shrinkage Mahalanobis structure. Keep it deployable:
weights must not use query truth.

Success criterion: learned weights beat unweighted and manual T-marker-weighted
baselines on held-out folds, Batch1-to-Batch2 transfer, and bootstrap stability.

### 6. Hierarchical or Refined OtherT Taxonomy

Split reference-side OtherT into interpretable sublabels such as double-negative
T, gamma-delta/MAIT-like, NKT-like, activated/cytotoxic-like, and Treg/Tfh-like
where marker evidence supports it. Transfer sublabels locally, then collapse
back to OtherT for package-level metrics.

Success criterion: collapsed OtherT precision and recall improve, and sublabels
show coherent marker profiles rather than unstable clustering artifacts.

### 7. Repeat-Frequency Estimator Tuning

Treat soft abundance as a compositional estimator. Tune `k`, bandwidth,
normalization, and optional OOD gating using true repeat samples. Report
bootstrap CIs, paired deltas, sign tests, and split-half ceiling.

Success criterion: median metacluster CLR weighted-r improves over 0.9566 or the
current soft estimator is confirmed as the locked default.

### 8. Runtime/API Hardening

Add chunked projection to query/reference constructors, preserve `label_prob` in
`somalign_query_from_som()`, reuse SOM assignments when available, stream
corrected projection, add compact storage modes, and expose `BPPARAM`.

Success criterion: large query/reference workflows reduce peak memory without
changing label, frequency, or correction outputs under equivalence tests.

## Recommended Immediate Sequence

1. Add a reproducible validation matrix script that reruns repeat frequencies,
   label transfer, calibration, and runtime smoke checks with structured outputs.
2. Validate high-margin OtherT rescue at 120k scale.
3. Build and test a strict biological T/NK router plus local specialist
   projection.
4. Add package-level calibration and soft-posterior APIs after the rescue/router
   experiment defines the required probability surfaces.
5. Address runtime correctness bugs and memory controls in a separate package
   patch, with equivalence tests.

## Release Claim Rules

Every performance statement in README, NEWS, vignettes, or manuscript text
should map to:

- dataset and split,
- command or script,
- package version and git SHA,
- source artifact hashes,
- primary and secondary metrics,
- raw prediction table,
- confidence interval or paired test where applicable,
- pass/fail threshold.

Until that exists, claims should be phrased as pilot validation results rather
than broad package guarantees.

## Bottom Line

SOMalign's core design is sound. The next large performance gain is most likely
to come from calibrated decision logic and local T/NK biology, while repeat
abundance should move toward soft posterior frequency estimation. Runtime work
should focus on projection memory and object retention, not replacing the OT
solver.
