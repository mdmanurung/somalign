# Do anchor (repeat) samples improve label transfer?

Experiment: `anchor_benefit_experiment.R` (120k stratified pilot cells, 6 lineages,
`gate_lineage2` ground truth, 15×15 SOM). We sweep the anchor cost-bonus strength
`rho_anchor` (0 = no anchors) and score transferred labels against the known label on
held-out (non-anchor) cells with `somalign_anchor_benefit()`.

## Headline

**Anchors can substantially rescue label transfer under severe, population-specific batch
effects — but only at `rho_anchor` two to three orders of magnitude above the default of 1,
and not at all when the batch effect is mild. At the production default `rho_anchor = 1` the
label-transfer benefit is nil.**

## Condition A — real acquisition-batch effect (pilot Batch1 → Batch2)

Plain label transfer already scores MCC 0.977 / accuracy 0.98 at 95% coverage; anchors add
essentially nothing (best lift +0.0002). Real same-study acquisition batches carry a mild,
largely global batch effect that optimal transport handles on its own. This mirrors the real
BMV result (label transfer works well without anchors) and explains why the production BMV fit,
run at `rho_anchor = 1`, saw no label benefit from its 10,000 anchors.

## Condition B — controlled severity sweep (per-lineage shift toward nearest neighbour)

Each lineage is moved toward its nearest neighbouring lineage by `delta` (δ = 1 ≈ full
collision), so plain OT progressively mis-maps populations. Anchor lift on MCC:

| δ (severity) | baseline (no anchors) | `rho_anchor = 1` | best | best `rho_anchor` | lift |
|---|---|---|---|---|---|
| 0.0 | 0.977 | 0.977 | 0.982 | 500 | +0.005 |
| 0.3 | 0.956 | 0.956 | 0.956 | 0 | 0.000 |
| 0.6 | 0.489 | 0.489 | 0.889 | 500 | **+0.400** |
| 0.9 | −0.044 | −0.044 | 0.537 | 500 | **+0.582** |
| 1.2 | −0.123 | −0.123 | 0.358 | 500 | **+0.481** |
| 1.5 | −0.225 | −0.225 | 0.389 | 500 | **+0.614** |

Two facts stand out. First, `rho_anchor = 1` is **identical to the no-anchor baseline at every
severity** — the default anchor strength is inert for label transfer. Second, large
`rho_anchor` (≈500) turns a collapsed transfer (MCC ≤ 0 at δ ≥ 0.9, i.e. worse than random)
into a partial recovery (MCC 0.36–0.54), because the cost bonus must be large enough to
overcome the cost gap that favours the wrong, collided mapping.

## Takeaways

- Anchors are worth using only when the batch effect is **population-specific and severe enough
  that plain OT mis-maps lineages**; they do nothing for mild or global shifts.
- The default `rho_anchor = 1` is far too weak to matter for labels. If anchors are intended to
  aid label transfer under a hard batch effect, `rho_anchor` should be tuned much higher
  (orders of magnitude), e.g. via a labelled held-out batch and `somalign_anchor_benefit()`.
- Even at best, anchors only partially recover severe collisions (MCC ~0.4–0.5, ~70% coverage),
  so they are a mitigation, not a fix, for badly confounded batches.

Caveat: the real BMV query is unlabelled by design, so anchor benefit cannot be measured on it
directly; this experiment uses the labelled pilot data (real batches + a controlled severity
sweep) as the closest measurable proxy.
