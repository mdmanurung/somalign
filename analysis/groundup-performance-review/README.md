# Ground-Up Performance Review

Date: 2026-07-19

This folder captures a multi-team review of SOMalign performance, starting from
the package's mathematical contract and ending in a ranked improvement program.
It is separate from `analysis/other-t-improvement/`, which contains the current
empirical BMV/OtherT experiment outputs.

## Inputs

- Core package code: `R/fit.R`, `R/ot.R`, `R/soft.R`, `R/results.R`,
  `R/query.R`, `R/reference.R`, `R/diagnostics.R`, `R/anchored.R`.
- Existing analysis folders:
  - `analysis/other-t-improvement/`
  - `analysis/label-validation/`
  - `analysis/anchor-benefit/`
  - `analysis/plans/2026-07-16-ideas-1-9/`
- Package readiness artifacts:
  - `somalign.Rcheck/00check.log`
  - `benchmarks/RESULTS.md`
  - `NEWS.md`

## Multi-Team Roles

- Team A, mathematical/statistical formulation: OT objective, uncertainty,
  calibration, soft posteriors, and metric learning.
- Team B, biological label transfer: OtherT/T/NK biology, routers, marker gates,
  label taxonomy, and deployable versus oracle boundaries.
- Team C, validation and publication readiness: repeat evidence, held-out label
  transfer, calibration, benchmarks, and release gates.
- Team D, runtime/API engineering: projection memory, SOM reuse, compact objects,
  chunking, solver and diagnostic runtime.

## Current Consensus

The dominant accuracy bottleneck is not the Sinkhorn solver. It is decision
control around a compressed SOM-node transfer: global `epsilon` trades OtherT
precision against recall, mixed T/NK nodes dilute node-level label probabilities,
and accepted-cell metrics can hide abstention costs.

The dominant reproducibility win is soft abundance. In the SOMalign example
notebook repeat metadata, 20 pilot-vs-BMV repeat pairs have median metacluster
CLR weighted-r of 0.9322 under hard nearest-node frequencies and 0.9566 under
soft kNN frequencies.

The dominant runtime bottleneck is nearest-code projection and retained dense
cell-level matrices, not OT over SOM nodes.

## Files

- `00_groundup_multiteam_synthesis.md` - final roundtable synthesis and ranked
  work program.
- `01_validation_matrix.md` - concrete experiments, success criteria, and
  release gates for the next improvement round.
