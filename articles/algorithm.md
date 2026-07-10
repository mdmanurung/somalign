# Algorithm: step-by-step

`somalign` aligns a query SOM to a fixed reference SOM using
codebook-level KL-unbalanced entropic optimal transport. The schematic
below traces each stage from raw data to the final per-sample output
columns.

    ╔══════════════════════════════════════════════════════════════════════════════════╗
    ║                        somalign · Algorithm Schematic                            ║
    ╚══════════════════════════════════════════════════════════════════════════════════╝

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 0 — TRAINING BOTH SOMs

      OLD data  (reference cohort)           NEW data  (query cohort)
      ┌──────────────────────────┐           ┌──────────────────────────┐
      │ T · N T · N T · M T · M │           │ · · · · · · · · · · · · │
      │ · T · · N N · M · · M · │           │ · · · · · · · · · · · · │
      │ T · · T · N · M M · · M │           │ · · · · · · · · · · · · │
      └──────────────────────────┘           └──────────────────────────┘
        labeled: T=T-cell  N=NK  M=Mono        unlabeled cells

        │ (1) z-score: compute μ, σ             │ (1) transform with REFERENCE μ, σ
        │ (2) train kohonen SOM                 │     (never new-data z-scores)
        ▼                                       │ (2) train kohonen SOM
                                                ▼
      REFERENCE SOM  C_ref  [K × p]          QUERY SOM  C_query  [M × p]
      ┌──────────────────────────────┐        ┌─────────────────────┐
      │ ◉TT  ◉TN  ◉NN  ◉NN  ◉NU    │        │  ◉q₁  ◉q₂  ◉q₃     │
      │ ◉TT  ◉TN  ◉NN  ◉NM  ◉MM    │        │  ◉q₄  ◉q₅  ◉q₆     │
      │ ◉MM  ◉MM  ◉MU  ◉UU  ◉UU    │        │  ◉q₇  ◉q₈  ◉q₉     │
      └──────────────────────────────┘        └─────────────────────┘
      K nodes. Each node stores:              M nodes. Each node stores:
      · centroid C_ref[k]  [1×p]             · centroid C_query[i]  [1×p]
      · label distribution                   · mass a[i] = fraction of new
      · mass b[j]                              samples mapped to this node
      · distance quantile threshold
      Labels show majority composition
      (UU = unlabeled/mixed)

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 1 — DIRECT PROJECTION  ★ PRIMARY RESULT  (no OT needed)

      For every new sample x_s, find the nearest reference SOM node in feature space:

      Feature 2 ▲
                │    ◉ r₁       ◉ r₂                   ◉ r₅
                │
                │                        × x_s
                │                   d₃↗     ↘d₄
                │                  ↗           ↘
                │    ◉ r₃        ◉               ◉ r₄
                │
                └──────────────────────────────────────▶ Feature 1

        k* = argmin_k  ‖ x_s − C_ref[k] ‖        d₃ < d₄  ⟹  k* = r₃

      Output per sample:
        old_som_unit            = k*
        old_som_distance        = ‖ x_s − C_ref[k*] ‖
        outside_reference_dist  = distance > quantile threshold?   ← novelty flag
        final_status            = "inside_reference" | "outside_reference" | "unknown"
        old_som_label           = majority label at k*
        old_som_label_confidence = fraction of k*'s training mass with that label

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 2 — OT ON CODEBOOKS  (KL-unbalanced Sinkhorn; runs once, not per sample)

      ┌─ 2a  Cost matrix D [M × K] ──────────────────────────────────────────────┐
      │                                                                            │
      │                  r₁    r₂    r₃    r₄    r₅                              │
      │           q₁  [ ░░    ████  ████  ████  ████ ]  ← q₁ is close to r₁     │
      │           q₂  [ ████  ████  ░░    ▒▒    ████ ]  ← q₂ is close to r₃     │
      │           q₃  [ ████  ████  ████  ████  ░░   ]  ← q₃ is close to r₅     │
      │                                                                            │
      │  D[i,j] = ‖ C_query[i] − C_ref[j] ‖²                                    │
      │  ░ low (cheap to transport)   ▒ medium   ████ high (expensive)            │
      └────────────────────────────────────────────────────────────────────────────┘

      ┌─ 2b  Solve for transport plan P [M × K] ─────────────────────────────────┐
      │                                                                            │
      │  min_{P≥0}  ΣΣ P[i,j] · D[i,j]       move mass along cheap paths         │
      │           + ε · KL(P ‖ a⊗b)          smooth / regularise the plan        │
      │           + ρ_q · KL(P1 ‖ a)         query may discard unmatched mass     │
      │           + ρ_r · KL(Pᵀ1 ‖ b)        ref may discard unmatched mass      │
      │                                                                            │
      │  KL-unbalanced: novel query nodes are not forced to match any ref node.   │
      │  Solved by alternating Sinkhorn scaling (pure-R or Python POT backend).   │
      └────────────────────────────────────────────────────────────────────────────┘

      ┌─ 2c  Resulting flows  (line weight ∝ P[i,j]) ────────────────────────────┐
      │                                                                            │
      │  Query nodes          mass flows              Reference nodes              │
      │                                                                            │
      │    ●─q₁  ═════════════════════════════════▶  r₁─◉  [T-cell 92%]          │
      │          ═══════════════════════════════▶     r₂─◉  [T-cell 78%]          │
      │                                                                            │
      │    ●─q₂          ══════════════════════▶     r₃─◉  [NK    88%]            │
      │                   ═════════════════════▶      r₄─◉  [NK    71%]            │
      │                                                                            │
      │    ●─q₃  · · · · · · · · · · · · · · ▶       r₅─◉  [Mono  95%]           │
      │           (tiny mass — novel node)                                         │
      │                                                                            │
      │  match_fraction[q₁] = Σⱼ P[1,j]/a[1] ≈ 0.95  → label transfer accepted  │
      │  match_fraction[q₃] = Σⱼ P[3,j]/a[3] ≈ 0.04  → label transfer rejected  │
      └────────────────────────────────────────────────────────────────────────────┘

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 3 — CORRECTION VECTORS  (per query node)

      For query node i, compute a weighted mean pull toward matched reference nodes:

             Δᵢ  =  Σⱼ  P[i,j] · (C_ref[j] − C_query[i])
                    ──────────────────────────────────────
                                Σⱼ P[i,j]

      Cartoon (2D feature space, node i = q₁):

      Feature 2 ▲
                │          ◉ r₁  ◀── P[1,r₁] large
                │         ↗
                │        ↗
                │   ◉ r₂  ◀─── P[1,r₂] smaller        pulls weighted by transport
                │      ↘ ↗                              mass sent to each ref node
                │       ●  q₁  ──────── Δ₁ ──────▶  ×  (shifted centroid)
                │
                └──────────────────────────────────────▶ Feature 1

      ‖Δᵢ‖ ≈ 0   node already aligns with its reference counterpart
      ‖Δᵢ‖ large  systematic displacement; batch effect or population shift

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 4 — CORRECTED PROJECTION  (auxiliary)

      Every sample in query node i is shifted by the same Δᵢ, then re-projected:

      Feature 2 ▲
                │      ◉ r_A              ◉ r_B
                │
                │            ★ x_s + Δᵢ  ════════════════▶ r_B  corrected_som_unit
                │           ↗
                │      Δᵢ ↗   (node-level shift; same for all samples in node i)
                │     ↗
                │    × x_s  ─────────────────────────────▶ r_A  old_som_unit
                │
                └──────────────────────────────────────────▶ Feature 1

        corrected_som_unit = argmin_k  ‖ (x_s + Δᵢ) − C_ref[k] ‖
        correction_norm    = ‖Δᵢ‖

      Sample reassigns r_A → r_B after OT correction.
      Both are reported; r_A (old_som_unit) is the primary, conservative result.

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      STAGE 5 — LABEL TRANSFER  (auxiliary)

      For each query node i, read row P[i, ·] to find the dominant reference node:

        dominant_j = argmax_j P[i,j]

      Bar chart of relative mass across reference nodes (node i = q₁):

        r₁  [T-cell 92%]  ████████████████████████  ← dominant (most mass)
        r₂  [T-cell 78%]  ████████████░░░░░░░░░░░░
        r₃  [NK    88%]   █████░░░░░░░░░░░░░░░░░░░
        r₄  [NK    71%]   ███░░░░░░░░░░░░░░░░░░░░░
        r₅  [Mono  95%]   ░░░░░░░░░░░░░░░░░░░░░░░░

      transferred_label         = label of r₁ = "T-cell"
      transferred_label_conf    = label purity at r₁ = 0.92

      Acceptance gate:
      ┌──────────────────────────────────────────────────────────────────────┐
      │                                                                      │
      │   match_fraction[i] ≥ 0.05      AND      confidence ≥ 0.60          │
      │                                                                      │
      │   PASS ─▶  transferred_label_accepted = TRUE                        │
      │   FAIL ─▶  transferred_label_accepted = FALSE                       │
      │            (node is novel or dominant ref node is too mixed)        │
      └──────────────────────────────────────────────────────────────────────┘

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

      OUTPUT COLUMNS  (one row per sample)

      ┌──────────────────────────────────┬─────────┬────────────────────────────────┐
      │ Column                           │ Stage   │ Purpose                        │
      ├──────────────────────────────────┼─────────┼────────────────────────────────┤
      │ old_som_unit                     │ 1       │ Primary node assignment         │
      │ old_som_distance                 │ 1       │ Distance to that node           │
      │ outside_reference_distance       │ 1       │ Novelty flag                   │
      │ final_status                     │ 1       │ inside / outside / unknown     │
      │ old_som_label                    │ 1       │ Primary label                  │
      │ old_som_label_confidence         │ 1       │ Label purity at node            │
      ├──────────────────────────────────┼─────────┼────────────────────────────────┤
      │ corrected_som_unit               │ 4       │ OT-corrected assignment        │
      │ correction_norm                  │ 3 / 4   │ Shift magnitude ‖Δᵢ‖           │
      ├──────────────────────────────────┼─────────┼────────────────────────────────┤
      │ transferred_label                │ 5       │ OT-derived label               │
      │ transferred_label_confidence     │ 5       │ Purity at dominant ref node    │
      │ transferred_label_accepted       │ 5       │ Gate result (TRUE / FALSE)      │
      └──────────────────────────────────┴─────────┴────────────────────────────────┘

      COMPLEXITY NOTE
      OT input is M × K codebook nodes, not n_new × n_ref individual samples.
      OT solve: O(M · K · iterations) — stays small even at n_new = 10⁶.
      Per-sample cost: O(n_new · K) nearest-node search, chunked in somalign_fit().
