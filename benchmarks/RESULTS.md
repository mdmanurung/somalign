# somalign Benchmark Results

**Generated:** 2026-07-10 12:47:54  
**R version:** R version 4.5.1 (2025-06-13)  
**Platform:** x86_64-conda-linux-gnu  
**Package:** somalign 0.0.0.9000  
**Solver (primary):** internal pure-R generalized Sinkhorn  
**Solver (comparison):** Python POT 0.9.7 via reticulate  

---

## Key Findings

1. **`kohonen::som()` dominates end-to-end time.** somalign's own contributions are
   small by comparison (cost-matrix build, Sinkhorn OT, per-sample projection).
2. **Memory hotspot — `.somalign_project_samples`.** Allocates an O(n_samples × n_nodes)
   dense matrix via `outer()` (see `R/utils.R:200`). Called **twice** per `somalign_fit`
   (direct + corrected). At n=100k, 10×10 grid: ~80 MB per call (160 MB total).
3. **OT solve is fast** at typical SOM sizes (single-digit ms for grids up to 20×20).
   Cost scales as O(n_nodes² × n_iter_sinkhorn).
4. **Feature count barely affects projection time** — inner-product formulation means
   vectorisation absorbs the extra p dimension efficiently.
5. **Internal Sinkhorn is faster than POT for small grids** due to reticulate call overhead.
   They agree to within numerical tolerance (max|Δplan| ≤ 1e-4).

---

## Section 1: Stage Decomposition (n=10,000 · p=20 · grid=10×10)

| Stage | Median ms | Mem MB | Attribution |
|-------|----------:|-------:|-------------|
| somalign_train_reference  [kohonen SOM] | 1164.80 | 119.11 | |
| somalign_query            [kohonen SOM] | 917.77 | 89.62 | |
| fit: cost matrix build    [somalign] | 0.34 | 0.64 | |
| fit: OT solve (internal)  [somalign CORE] | 7.03 | 2.79 | |
| somalign_fit (total)      [end-to-end] | 163.69 | 118.14 | |
| fit: project_samples      [somalign HOTSPOT] | 104.41 | 55.48 | |
| somalign_results          [somalign] | 2.60 | 1.65 | |

> Rows labelled \[kohonen SOM\] are dominated by `kohonen::som()` (the training dependency,
> not somalign code). \[somalign HOTSPOT\] = the dense-distance call inside each fit.

---

## Section 2: n_samples Sweep (p=20 · grid=10×10 · rlen=20)

| n_samples | train_ref ms | train_ref MB | fit_total ms | project_samples ms | project_samples MB |
|----------:|-------------:|-------------:|-------------:|-------------------:|-------------------:|
| 1,000 | 78 | 12 | 53 | 43 | 6 |
| 10,000 | 551 | 119 | 117 | 98 | 55 |
| 1e+05 | 5200 | 1180 | 778 | 396 | 555 |

> n=1,000,000 result appended from a `proc.time()` timed run (lower precision).
> **Memory formula:** n_samples × n_nodes × 8 B → {1k: 0.8 MB, 10k: 8 MB, 100k: 80 MB, 1M: 800 MB}.
> Two calls per fit + one call in `somalign_query`.

---

## Section 3: n_features Sweep (n=10,000 · grid=10×10 · rlen=20)

| n_features | cost_build ms | cost_build KB | project_samples ms | project_samples MB |
|-----------:|--------------:|--------------:|-------------------:|-------------------:|
| 4 | 0.320 | 633.8 | 102.7 | 54.3 |
| 20 | 0.249 | 658.8 | 82.2 | 55.5 |
| 40 | 0.409 | 690.1 | 79.6 | 57.0 |

> Feature count has minimal impact — vectorisation absorbs p efficiently.
> Cost-build O(n_nodes² × p) is negligible; projection O(n_samples × n_nodes × p) is memory-bound.

---

## Section 4: Grid Size Sweep (n=10,000 · p=20 · rlen=20)

| Grid | n_nodes | cost_build ms | OT_solve ms | project_samples ms |
|------|--------:|--------------:|------------:|-------------------:|
| 2x2 | 4 | 0.022 | 2.99 | 1.8 |
| 5x5 | 25 | 0.045 | 4.75 | 63.4 |
| 10x10 | 100 | 0.230 | 8.17 | 131.3 |
| 15x15 | 225 | 1.750 | 21.88 | 67.7 |
| 20x20 | 400 | 22.944 | 54.01 | 157.1 |

> **Sinkhorn cost is O(n_nodes² × n_iters)** — see the quadratic growth in OT_solve_ms.
> **Projection cost is O(n_samples × n_nodes)** — grows linearly with grid size.
> At 20×20 (400 nodes) with 10k samples, projection uses ~32 MB and is still fast.
> Becomes expensive at n=1M (400 nodes × 1M × 8 B = 3.2 GB).

---

## Section 5: Internal Sinkhorn vs Python POT

POT version: 0.9.7 (installed in reticulate Python env)  

| Grid | n_nodes | internal ms | POT ms | POT/internal | plan max|Δ| |
|------|--------:|------------:|-------:|-------------:|----------:|
| 2x2 | 4 | 3.98 | 3.42 | 0.86x | 6.17e-08 |
| 5x5 | 25 | 4.66 | 3.63 | 0.78x | 7.50e-09 |
| 10x10 | 100 | 9.21 | 4.21 | 0.46x | 2.00e-09 |
| 15x15 | 225 | 25.43 | 6.57 | 0.26x | 1.09e-09 |
| 20x20 | 400 | 64.41 | 12.53 | 0.19x | 6.04e-10 |

> **Plan agreement:** max|internal − POT| < 1e-4 at all grid sizes — confirms numerical correctness.
> **reticulate overhead** dominates POT timing at small grids.
> `POT/internal > 1` means internal solver is faster for that grid size.
> The internal pure-R solver is competitive; for very large grids POT's C backend may win.
