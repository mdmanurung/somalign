test_that("somalign_mnn_novelty: held-out-novel case — novel-population nodes flagged, shared-pop nodes clean", {
  # Reference has 2 populations (pop A near -4, pop B near +4) in 3D marker space.
  # Query has those 2 PLUS a genuinely novel population (near +20) far from both.
  # No batch shift on the shared populations so shared-pop nodes are clean at k=1.

  set.seed(7)
  p <- 3L
  n_per_pop <- 80L

  # Reference data: two compact clusters, well separated
  pop_A <- matrix(rnorm(n_per_pop * p, mean = -4, sd = 0.4), ncol = p)
  pop_B <- matrix(rnorm(n_per_pop * p, mean =  4, sd = 0.4), ncol = p)
  ref_data <- rbind(pop_A, pop_B)
  colnames(ref_data) <- paste0("M", seq_len(p))
  labels_ref <- c(rep("A", n_per_pop), rep("B", n_per_pop))

  # Query data: same two clusters (no shift) + novel population far away
  pop_A_q  <- matrix(rnorm(n_per_pop * p, mean = -4, sd = 0.4), ncol = p)
  pop_B_q  <- matrix(rnorm(n_per_pop * p, mean =  4, sd = 0.4), ncol = p)
  pop_C_q  <- matrix(rnorm(50L    * p, mean = 20, sd = 0.4), ncol = p)  # novel
  qry_data <- rbind(pop_A_q, pop_B_q, pop_C_q)
  colnames(qry_data) <- paste0("M", seq_len(p))
  labels_qry <- c(rep("A", n_per_pop), rep("B", n_per_pop), rep("C", 50L))

  g <- kohonen::somgrid(3L, 3L, "hexagonal")
  ref <- somalign_train_reference(ref_data, labels = labels_ref, grid = g, rlen = 30)
  qry <- somalign_query(qry_data, ref, grid = g, rlen = 30)
  fit <- somalign_fit(qry, ref)

  flags <- somalign_mnn_novelty(fit)

  # flags should be a named logical vector with one entry per query node
  expect_true(is.logical(flags))
  expect_equal(length(flags), nrow(fit$query$codebook))
  expect_equal(names(flags), paste0("q", seq_len(nrow(fit$query$codebook))))

  # Identify which query nodes are dominated by each population by tallying
  # cell-to-node assignments
  su <- fit$query$sample_unit
  node_label_tbl <- table(su, factor(labels_qry[seq_len(length(su))]))
  node_dom <- apply(node_label_tbl, 1, function(r) {
    if (sum(r) == 0) NA_character_ else colnames(node_label_tbl)[which.max(r)]
  })

  novel_nodes  <- as.integer(names(node_dom)[!is.na(node_dom) & node_dom == "C"])
  shared_nodes <- as.integer(names(node_dom)[!is.na(node_dom) & node_dom %in% c("A", "B")])

  # All novel-dominated nodes must be flagged
  expect_true(length(novel_nodes) > 0, label = "at least one novel-dominated node exists")
  expect_true(all(flags[novel_nodes]),
              label = "all novel-dominated nodes must be MNN-unmatched")

  # Shared-population nodes should mostly NOT be flagged (low false-fire rate)
  # We use a threshold rather than exact-zero because independently-trained SOMs
  # can have incidental non-reciprocity from unequal node allocation (pigeonhole).
  if (length(shared_nodes) > 0) {
    false_fire_rate_shared <- mean(flags[shared_nodes])
    expect_lt(false_fire_rate_shared, 0.4,
              label = "shared-pop nodes: false-fire rate < 40% at zero shift")
  }
})


test_that("somalign_mnn_novelty: batch-shift confound — false-fire rate rises with shift magnitude", {
  # Two shared populations, no novel population.
  # We measure three conditions, all on POPULATED NODES ONLY (nodes whose
  # dominant-population tally is non-empty), to exclude structurally-interstitial
  # nodes that inflate any rate regardless of shift.
  #
  #   (a) SANITY CHECK: identical-codebook baseline (reuse ref$som_ref as query
  #       SOM) → rate = 0 by mathematical identity (self is own nearest neighbour)
  #   (b) ZERO-SHIFT BASELINE: independently-trained query SOM, same data, no
  #       batch translation → measures the pure pigeonhole non-reciprocity baseline
  #       (~40–50% on populated nodes with 9 nodes / 2 populations).  This IS the
  #       confound E2 sees: no shift, but structural non-reciprocity exists.
  #   (c) LARGE PER-POPULATION SHIFT: one coherent translation vector per
  #       population applied uniformly to all its cells → rate rises further,
  #       documents that shift adds confound on top of the structural baseline.
  #
  # Key design choice: each population gets a SINGLE per-marker shift drawn from
  # N(0, shift_sd) applied to EVERY cell in that population.  This is a coherent
  # batch translation (not per-cell variance inflation) that keeps cluster
  # structure intact and exercises the confound the docstring warns about.

  set.seed(13)
  p          <- 4L
  n_per_pop  <- 100L
  shift_sd   <- 4     # coherent translation: 10 × within-cluster SD, clusters at ±5

  pop_A <- matrix(rnorm(n_per_pop * p, mean = -5, sd = 0.4), ncol = p)
  pop_B <- matrix(rnorm(n_per_pop * p, mean =  5, sd = 0.4), ncol = p)
  ref_data <- rbind(pop_A, pop_B)
  colnames(ref_data) <- paste0("M", seq_len(p))
  labels_ref <- c(rep("A", n_per_pop), rep("B", n_per_pop))

  g   <- kohonen::somgrid(3L, 3L, "hexagonal")
  ref <- somalign_train_reference(ref_data, labels = labels_ref, grid = g, rlen = 30)

  # Helper: compute false-fire rate on POPULATED nodes only via label tally
  populated_false_fire <- function(fit, flags, labels) {
    su <- fit$query$sample_unit
    node_label_tbl <- table(su, factor(labels))
    node_dom <- apply(node_label_tbl, 1, function(r) {
      if (sum(r) == 0) NA_character_ else colnames(node_label_tbl)[which.max(r)]
    })
    pop_nodes <- as.integer(names(node_dom)[!is.na(node_dom)])
    if (length(pop_nodes) == 0) return(NA_real_)
    mean(flags[pop_nodes])
  }

  # (a) Sanity: identical codebooks → reciprocity is trivially exact, rate = 0
  qry_sanity <- somalign_query(
    ref_data, ref,
    som_query = ref$som_ref, codebook_space = "reference_scaled"
  )
  fit_sanity <- somalign_fit(qry_sanity, ref)
  flags_sanity  <- somalign_mnn_novelty(fit_sanity)
  rate_sanity   <- mean(flags_sanity)   # no need for populated-only; all 0
  message(sprintf("MNN false-fire (identical codebook sanity): %.3f", rate_sanity))
  expect_equal(rate_sanity, 0,
               label = "identical-codebook sanity: every node is reciprocal → rate = 0")

  # (b) Zero-shift baseline: independently trained query SOM, same data, no translation
  set.seed(17)
  qry_zero <- somalign_query(ref_data, ref, grid = g, rlen = 30)
  fit_zero  <- somalign_fit(qry_zero, ref)
  flags_zero  <- somalign_mnn_novelty(fit_zero)
  rate_zero   <- populated_false_fire(fit_zero, flags_zero, labels_ref)
  message(sprintf("MNN false-fire (zero shift, populated nodes): %.3f  [pigeonhole baseline]", rate_zero))
  # Do not assert a floor; merely document.  The value (~0.4–0.5) is the
  # structural confound that E2 must account for even at SD = 0.
  expect_true(is.numeric(rate_zero) && !is.na(rate_zero),
              label = "zero-shift rate is a finite number")

  # (c) Large coherent per-population shift: one translation vector per pop
  set.seed(13)
  delta_A <- rnorm(p, 0, shift_sd)   # one vector for all pop-A cells
  delta_B <- rnorm(p, 0, shift_sd)   # one vector for all pop-B cells
  qry_data_large <- rbind(
    pop_A + matrix(delta_A, n_per_pop, p, byrow = TRUE),
    pop_B + matrix(delta_B, n_per_pop, p, byrow = TRUE)
  )
  colnames(qry_data_large) <- paste0("M", seq_len(p))
  qry_large <- somalign_query(qry_data_large, ref, grid = g, rlen = 30)
  fit_large <- suppressWarnings(somalign_fit(qry_large, ref))
  flags_large <- somalign_mnn_novelty(fit_large)
  rate_large  <- populated_false_fire(fit_large, flags_large, labels_ref)
  message(sprintf("MNN false-fire (coherent shift SD=%g, populated nodes): %.3f", shift_sd, rate_large))

  # A coherent population-level shift displaces whole clusters → substantially
  # elevated false-fire rate vs. the zero-shift baseline.
  expect_true(is.numeric(rate_large) && !is.na(rate_large),
              label = "large-shift rate is a finite number")
  expect_gt(rate_large, rate_zero,
            label = "coherent per-population batch shift raises false-fire rate above zero-shift baseline")
})


test_that("somalign_mnn_novelty: structural guards", {
  # Not a somalign_fit object
  expect_error(somalign_mnn_novelty(list()), "`fit` must be a `somalign_fit` object")
  expect_error(somalign_mnn_novelty("hello"), "`fit` must be a `somalign_fit` object")

  # NULL codebook guard: craft a minimal somalign_fit-like object
  fake_fit <- structure(
    list(
      query     = structure(list(codebook = NULL, sample_unit = 1L), class = "somalign_query"),
      reference = structure(list(codebook = matrix(1:4, 2, 2)), class = "somalign_reference")
    ),
    class = "somalign_fit"
  )
  expect_error(somalign_mnn_novelty(fake_fit), "fit\\$query\\$codebook.*NULL")
})


test_that("somalign_mnn_novelty_cells: broadcasts per-node flag to cells correctly", {
  set.seed(3)
  p <- 2L
  mat <- matrix(rnorm(60 * p), ncol = p, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  flags <- somalign_mnn_novelty(fit)
  cell_flags <- somalign_mnn_novelty_cells(fit, flags)

  # One logical per cell
  expect_equal(length(cell_flags), nrow(mat))
  expect_true(is.logical(cell_flags))

  # Manual broadcast: flags indexed by sample_unit
  expected <- as.logical(flags[fit$query$sample_unit])
  expect_identical(cell_flags, expected)
})


test_that("somalign_mnn_novelty_cells: structural guards", {
  expect_error(somalign_mnn_novelty_cells(list(), logical(4)), "`fit` must be a `somalign_fit` object")

  set.seed(5)
  p <- 2L
  mat <- matrix(rnorm(40 * p), ncol = p, dimnames = list(NULL, c("A", "B")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)

  # Wrong length
  expect_error(somalign_mnn_novelty_cells(fit, logical(1)), "logical vector of length equal to the number of query nodes")
  # Wrong type
  expect_error(somalign_mnn_novelty_cells(fit, rep(1L, nrow(fit$query$codebook))),
               "logical vector of length equal to the number of query nodes")
})
