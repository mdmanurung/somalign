# Statistical Physicist — 2 ideas

---

# Epsilon Phase-Transition Diagnostic and Critical-Epsilon Estimator

## Persona
**Statistical Physicist** — epsilon is temperature; the transport plan is a Gibbs
distribution whose "magnetisation" undergoes a continuous phase transition as
epsilon crosses a critical value.

## Motivation
A statistical physicist recognises that the UOT objective
`F(P) = <C, P> - epsilon * H(P) + rho_query * KL(a|P1) + rho_ref * KL(b|P^T 1)`
is a free-energy functional. The Gibbs plan `P_ij ∝ exp(-(C_ij - f_i - g_j) / epsilon)`
is the equilibrium state at temperature `epsilon`. In the zero-temperature limit
the plan concentrates on the optimal coupling; at infinite temperature it spreads
uniformly. The transition between a "localised" plan (each query node maps to one
reference node) and a "delocalised/collapsed" plan (all query mass pools) is
generically a continuous phase transition in epsilon — analogous to the
paramagnetic-to-ferromagnetic transition — with an order parameter that jumps
(or curves sharply) at a critical `epsilon_c`.

The current `somalign_sensitivity_grid` runs a full `.somalign_align_transport`
call for every (epsilon, rho) combination. This is expensive and gives no
mechanistic hint about *where* the interesting transition is. A physicist would
compute an order parameter along an epsilon sweep cheaply, locate the
susceptibility peak, and call that `epsilon_c` — the principled lower bound
below which the plan is transport-cost-dominated and above which it smears.

## Connection to Existing Code/Data

- `ot.R`: `.somalign_solve_internal` and `.somalign_solve_internal_log` already
  return the converged potentials `f`, `g` (log-domain) or `u`, `v` (primal).
  The log-partition function `log Z_i = logsumexp_j((f_i + g_j - C_ij) / eps)`
  is already computed row-by-row inside the log-domain Sinkhorn (the `lse_g`
  vector on line 123 of `ot.R`). Only the scalar `log Z = sum_i log Z_i` needs
  to be accumulated.
- `fit.R`: `.somalign_align_transport` already normalises the cost by
  `cost_scale` (median positive entry) so epsilon is already dimensionless.
- `diagnostics.R`: `somalign_sensitivity_grid` sweeps epsilon values and
  collects per-fit scalars into a data frame — the ideal framework to reuse.
- `fit$transport_plan` (the M×K matrix `P`): the plan is already stored, so
  the order parameter can be computed post-hoc from an existing fit without
  re-running Sinkhorn.

## Approach

1. **Define the order parameter.** For each query node `i`, the "localisation"
   is the effective number of reference nodes it uses:
   `S_i = exp(H(P_i))` where `H(P_i) = -sum_j p_ij * log(p_ij)` and `p_ij` is
   the row-normalised plan (already computed as `correspondence` in
   `.somalign_align_transport`). The system-level order parameter is
   `Phi(epsilon) = mean_i(S_i) / K`, the mean fractional usage of reference
   nodes. At `epsilon → 0`, `Phi → 1/K` (one node used); at `epsilon → inf`,
   `Phi → 1` (all nodes used equally). The "susceptibility" is
   `chi(epsilon) = d Phi / d epsilon` — its peak marks `epsilon_c`.

2. **Compute the free energy.** From the converged log-domain potentials expose
   `log Z = sum_i [epsilon * log a_i + f_i] + sum_j [epsilon * log b_j + g_j]`
   (the dual objective). This is the negative Helmholtz free energy times
   epsilon. Exposing `log Z` requires a one-line addition to
   `.somalign_solve_internal_log` to sum `lse_g` (already on the stack) and
   return it alongside `plan`. Add `free_energy` to the `diagnostics$solver`
   list in `.somalign_build_diagnostics`.

3. **Cheap epsilon sweep.** Add a new function
   `somalign_epsilon_sweep(query, reference, epsilon_grid, ...)` that runs
   `.somalign_solve_ot` across a log-spaced epsilon grid and records
   `(epsilon, Phi, log_Z, iterations, converged)` without the expensive
   per-cell projection step (skip `.somalign_project_pair`). Cost: one Sinkhorn
   per epsilon value.

4. **Critical-epsilon estimator.** From the sweep, estimate `epsilon_c` as the
   epsilon at which the numerical derivative `d Phi / d epsilon` (finite
   difference over the log-epsilon grid) is maximised. Report it alongside a
   plot of `Phi` vs `log(epsilon)` using the existing ggplot2 infrastructure
   (`R/plot.R`). Suggest `epsilon = 0.3 * epsilon_c` as the default working
   temperature (analogy: working below Curie temperature to stay in the ordered
   phase with headroom for numerical stability).

## Expected Improvement

- **Principled epsilon selection**: replaces the current manual
  `somalign_sensitivity_grid` with a data-driven recommendation. For the
  39.8M-cell BMV-to-pilot alignment the critical epsilon would have flagged
  whether the default 0.1 was already in the delocalised regime.
- **Free-energy diagnostic**: `log Z` is a single scalar that quantifies
  alignment quality independent of post-hoc projections. A large drop in
  `log Z` between fits signals a better-constrained transport.
- **Speed**: the sweep avoids per-cell projection entirely, so it runs ~10x
  faster than `somalign_sensitivity_grid` for the same epsilon range.

## Feasibility
- **Effort**: Medium
- **Fits current architecture**: Yes — `somalign_epsilon_sweep` is a new thin
  wrapper around `.somalign_solve_ot`; `free_energy` is a one-line addition to
  the log-domain solver.
- **Methods available**: Standard — entropy of discrete distributions, numerical
  differentiation, logsumexp (already implemented as `.somalign_logsumexp`).
- **Key risk**: The phase transition may be "washed out" (very gradual) for
  small codebooks (e.g. 4×4 SOM) where K is small and the plan is always
  moderately localised, making `epsilon_c` hard to locate. Validation on the
  real BMV dataset (large K, many nodes) should show a clear inflection.

---

# Simulated Annealing Sinkhorn: Epsilon Cooling Schedule as a Solver Option

## Persona
**Statistical Physicist** — the fixed-epsilon Sinkhorn iteration is isothermal
equilibration; annealing from high to low temperature should find the ground
state faster and with fewer local-minima artefacts than a single cold solve.

## Motivation
The current solvers (`.somalign_solve_internal`, `.somalign_solve_internal_log`)
fix epsilon throughout and iterate until potentials converge. This is isothermal.
A physicist knows that for rugged free-energy landscapes (here: high cost
heterogeneity, sparse masses, or label-guided cost penalties with huge entries),
isothermal solves at low epsilon get trapped in local minima or require many
more iterations.

Simulated annealing for OT means: start at a high epsilon (flat landscape, easy
to equilibrate), then slowly cool to the target epsilon, carrying warm potentials
forward as the initial condition for each subsequent temperature. Because the
Sinkhorn iterates are continuous in epsilon, the warm-start dramatically reduces
the number of iterations needed at each new temperature compared to cold-starting.
This is directly analogous to the "adiabatic" cooling schedule used in spin-glass
theory to avoid metastability.

The two-pass solver (`somalign_fit_two_pass`) hints at this philosophy — it uses
a high `epsilon_global` followed by a low `epsilon_local`. Annealing generalises
this to an arbitrary schedule with warm potentials carried between steps, which
is both faster and yields a better final plan.

## Connection to Existing Code/Data

- `ot.R`: `.somalign_solve_internal_log` already works with log-potentials
  `f` (M-vector) and `g` (K-vector) that can be passed as initial conditions.
  Currently they are initialised to zeros (lines 110-111 of `ot.R`). The only
  change needed is an `f_init` / `g_init` argument to warm-start from a
  previous temperature's converged potentials.
- `.somalign_sinkhorn_kernel` in the primal solver constructs `K = exp(-C/eps)`;
  annealing requires rebuilding K at each temperature, which is just one
  `exp(-cost / eps_new)` call — cheap.
- `fit.R`: `somalign_fit_two_pass` already demonstrates two-temperature
  warm-starting philosophically; the annealing solver is a principled
  generalisation with an arbitrary cooling schedule.
- The `label_guided` cost mask in `.somalign_prepare_cost` (line 199: `penalty
  = max(cost_normalized) * 1e4`) creates a very rugged landscape — the exact
  regime where annealing beats isothermal.
- `fit.R`: `somalign_fit` and `.somalign_solve_ot` accept `solver` as a
  character argument, making it easy to add `"annealing"` as a new option
  alongside `"internal"` and `"log_domain"`.

## Approach

1. **Warm-start interface.** Modify `.somalign_solve_internal_log` to accept
   optional `f_init` and `g_init` arguments (default `NULL`, falling back to
   zero initialisation). Add a parameter check; no other changes to the
   convergence loop.

2. **Cooling schedule.** Implement a helper
   `.somalign_cooling_schedule(epsilon_start, epsilon_target, n_steps, type)`
   returning a decreasing sequence. Default: geometric schedule
   `eps_k = epsilon_start * (epsilon_target / epsilon_start)^(k / n_steps)`,
   `k = 0, ..., n_steps`. Default `n_steps = 10`, `epsilon_start = 10 *
   epsilon_target`. Both parameters are user-tunable.

3. **Annealing solver.** Add `.somalign_solve_annealing(cost, a, b, epsilon,
   rho_query, rho_ref, max_iter, tol, n_steps, epsilon_start)` that: (a) builds
   the schedule; (b) runs `.somalign_solve_internal_log` at `eps_k` for
   `max_iter / n_steps` iterations; (c) passes converged `f`, `g` as warm start
   to the next step; (d) at the last step, runs to full convergence tolerance.
   Total iteration budget is identical to the single-solve budget, split across
   temperatures.

4. **Expose as `solver = "annealing"`.** Add `"annealing"` to the `solver`
   `match.arg` in `somalign_fit` (`fit.R` line 97) and dispatch to
   `.somalign_solve_annealing` in `.somalign_solve_ot` (`ot.R` lines 14-18).
   Expose `n_steps` and `epsilon_start` as optional arguments to `somalign_fit`
   (with defaults so the interface is backward-compatible). Record the cooling
   schedule in `diagnostics$solver` so `somalign_diagnostics` surfaces it.

## Expected Improvement

- **Convergence speed for hard problems**: for label-guided fits (`label_guided =
  TRUE`) with the `1e4 * max_cost` penalty, annealing is expected to converge in
  fewer total Sinkhorn iterations because the warm potentials at each temperature
  are already close to the new equilibrium. Benchmark against
  `somalign_sensitivity_grid` on the BMV/pilot data.
- **Better ground state for rugged cost landscapes**: the annealing plan should
  have lower primal objective value than the cold-start plan at the same final
  epsilon, which means tighter barycentric node shifts in `.somalign_node_shifts`
  and less shrinkage artefact in the corrected projection.
- **Graceful handling of small epsilon without log-domain solver**: annealing
  reaches small epsilon from above, keeping the primal kernel `K = exp(-C/eps)`
  numerically safe at intermediate temperatures, so the standard solver can be
  used deeper into the small-epsilon regime before the log-domain solver is
  needed.

## Feasibility
- **Effort**: Low
- **Fits current architecture**: Yes — the warm-start interface is a two-line
  change to `.somalign_solve_internal_log`; the annealing driver and schedule
  helper are new ~40-line functions; the dispatch hook is a `switch` case in
  `.somalign_solve_ot`.
- **Methods available**: Standard — geometric cooling schedules and warm-started
  Sinkhorn are textbook methods in computational OT (Peyré & Cuturi 2019,
  Chapter 4).
- **Key risk**: With a fixed total iteration budget split across `n_steps`
  temperatures, each temperature may not reach its local equilibrium before
  cooling, which can leave the final solve starting from a suboptimal warm point.
  Mitigation: use a per-temperature inner convergence check and only cool when
  `delta < tol_anneal` (a looser tolerance, e.g. `10 * tol`), falling back to
  the full `max_iter` if needed.
