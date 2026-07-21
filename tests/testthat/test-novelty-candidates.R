# Tests for somalign_novelty_candidates()
#
# Design principles:
# - Novel populations are placed far from the reference (large displacement, ~8
#   units) so their k-NN mean distance scores cleanly top the distribution and
#   dominate the high-novelty tail.
# - tail_quantile and min_cluster are tuned per test: the tail must contain
#   >= min_cluster novel cells per group.  With novel_center at distance ~8 from
#   origin and in-distribution cells near origin, a sufficiently low
#   tail_quantile (0.65-0.80) captures all novel cells and rejects in-dist.
# - Artifact populations are placed at DIFFERENT reference-scaled coordinates in
#   each group; because they do not recur at consistent coordinates, they are not
#   linked across groups and are therefore not minted.
# - withr::local_seed() is used to make k-means deterministic.
#
# All tests work in 2D reference-scaled space (center=0, scale=1).

# ---------------------------------------------------------------------------
# Shared fixture builder
# ---------------------------------------------------------------------------

# Build a synthetic somalign_fit with:
#   - reference: 4 nodes on a 2x2 grid near origin
#   - query:     n_ref_per_group * n_groups in-distribution cells +
#                optional novel population at `novel_center` present in all groups
#                optional per-group artifact at different centers
make_novelty_fixture <- function(
    seed              = 42L,
    n_groups          = 3L,
    n_ref_per_group   = 200L,
    novel_center      = c(8, 0),
    n_novel_per_group = 80L,
    artifact_centers  = NULL,  # list of length n_groups; each is c(x,y)
    n_artifact_per_group = 0L
) {
  withr::local_seed(seed)

  # Reference: 4 nodes at corners of a small square near origin
  ref_codebook <- matrix(
    c(-0.5, -0.5,
       0.5, -0.5,
      -0.5,  0.5,
       0.5,  0.5),
    nrow = 4, byrow = TRUE
  )
  colnames(ref_codebook) <- c("x", "y")

  ref <- somalign_reference_from_nodes(
    codebook    = ref_codebook,
    features    = c("x", "y"),
    center      = c(x = 0, y = 0),
    scale       = c(x = 1, y = 1),
    node_masses = rep(0.25, 4),
    label_prob  = matrix(
      c(1, 0,
        1, 0,
        0, 1,
        0, 1),
      nrow = 4, byrow = TRUE,
      dimnames = list(NULL, c("A", "B"))
    ),
    distance_quantiles = matrix(
      rep(c(0.5, 1.0, 1.5, 2.0), 4),
      nrow = 4, byrow = TRUE,
      dimnames = list(NULL, c("50%", "90%", "95%", "99%"))
    )
  )

  groups_vec <- character(0)
  query_rows <- matrix(numeric(0), ncol = 2)
  colnames(query_rows) <- c("x", "y")

  for (g in seq_len(n_groups)) {
    grp_name <- paste0("batch", g)

    # In-distribution cells near origin
    in_dist <- matrix(
      c(rnorm(n_ref_per_group, 0, 0.3),
        rnorm(n_ref_per_group, 0, 0.3)),
      ncol = 2
    )
    colnames(in_dist) <- c("x", "y")
    query_rows <- rbind(query_rows, in_dist)
    groups_vec <- c(groups_vec, rep(grp_name, n_ref_per_group))

    # Novel cells: same location in all groups
    if (n_novel_per_group > 0) {
      nov <- matrix(
        c(rnorm(n_novel_per_group, novel_center[1], 0.1),
          rnorm(n_novel_per_group, novel_center[2], 0.1)),
        ncol = 2
      )
      colnames(nov) <- c("x", "y")
      query_rows <- rbind(query_rows, nov)
      groups_vec <- c(groups_vec, rep(grp_name, n_novel_per_group))
    }

    # Per-group artifact: different location each group
    if (!is.null(artifact_centers) && n_artifact_per_group > 0) {
      ac <- artifact_centers[[g]]
      art <- matrix(
        c(rnorm(n_artifact_per_group, ac[1], 0.1),
          rnorm(n_artifact_per_group, ac[2], 0.1)),
        ncol = 2
      )
      colnames(art) <- c("x", "y")
      query_rows <- rbind(query_rows, art)
      groups_vec <- c(groups_vec, rep(grp_name, n_artifact_per_group))
    }
  }

  # SOM for query: reference nodes + novel node
  query_codebook <- rbind(
    ref_codebook,
    matrix(c(novel_center[1], novel_center[2]), nrow = 1,
           dimnames = list(NULL, c("x", "y")))
  )
  query_som <- make_som(query_codebook)

  qry <- somalign_query(query_rows, ref, som_query = query_som)
  fit <- somalign_fit(qry, ref, solver = "internal")

  list(fit = fit, groups = groups_vec, novel_center = novel_center,
       ref = ref, query_rows = query_rows,
       n_groups = n_groups,
       n_ref_per_group = n_ref_per_group,
       n_novel_per_group = n_novel_per_group)
}


# ---------------------------------------------------------------------------
# Test 1: Reproducible novel population IS minted
# ---------------------------------------------------------------------------
test_that("reproducible novel population is minted and graftable", {
  withr::local_seed(1L)

  # 3 groups x 200 in-dist + 100 novel each = 300 per group = 900 total.
  # Novel cells score ~8; in-dist ~0.7.
  # tail_quantile=0.65 -> top 35% = 315 cells.
  # All 300 novel cells score ~8 -> top 300; plus 15 in-dist cells.
  # Per group: ~100 novel + ~5 in-dist = ~105 tail cells.
  # min_cluster=60 -> 105 %/% 60 = 1; k_g=1 -> one cluster of 105 cells >= 60. Minted.
  fix <- make_novelty_fixture(
    seed              = 101L,
    n_groups          = 3L,
    n_ref_per_group   = 200L,
    novel_center      = c(8, 0),
    n_novel_per_group = 100L
  )

  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile = 0.65,
    min_cluster   = 60L,
    min_batches   = 2L,
    tol_factor    = 2.0
  )

  expect_s3_class(cand, "somalign_novelty_candidates")

  # At least one candidate minted
  expect_gte(nrow(cand$prototypes), 1L)

  # Prototype lands near the true novel centroid (within 1.5 units)
  novel_proto <- cand$prototypes[1, ]
  dist_to_novel <- sqrt(sum((novel_proto - fix$novel_center)^2))
  expect_lt(dist_to_novel, 1.5)

  # n_groups_support >= 2
  expect_gte(cand$n_groups_support[1], 2L)

  # Full-length score and tail
  n_cells <- nrow(fix$fit$query$scaled_data)
  expect_length(cand$score, n_cells)
  expect_length(cand$tail,  n_cells)

  # Novel cells score much higher than in-dist cells
  novel_idx   <- which(abs(fix$query_rows[, "x"] - 8) < 1.0)
  in_dist_idx <- which(abs(fix$query_rows[, "x"]) < 0.5 &
                         abs(fix$query_rows[, "y"]) < 0.5)
  expect_gt(median(cand$score[novel_idx]),
            median(cand$score[in_dist_idx]))

  # Graft: extend reference with prototypes
  extended <- somalign_extend_reference(
    fix$ref,
    cand$prototypes,
    new_labels = paste0("novel_", seq_len(nrow(cand$prototypes)))
  )
  expect_s3_class(extended, "somalign_reference")
  expect_gte(nrow(extended$codebook), nrow(fix$ref$codebook) + 1L)

  # Refit and check novel label assignment.
  # Row order: somalign_results() returns one row per query cell in the same
  # order as fix$query_rows (query$sample_id is a sequential index).
  new_query <- somalign_query(
    fix$query_rows, extended,
    som_query = make_som(rbind(fix$ref$codebook, cand$prototypes))
  )
  new_fit <- somalign_fit(new_query, extended, solver = "internal")
  results  <- somalign_results(new_fit)

  # Verify row alignment: one result row per query cell.
  expect_equal(nrow(results), nrow(fix$query_rows))

  # Identify the known novel cells by their coordinate (|x - 8| < 1).
  novel_idx <- which(abs(fix$query_rows[, "x"] - 8) < 1.0)
  expect_gt(length(novel_idx), 0L)

  # The direct projection label (old_som_label) is independent of the
  # acceptance gate: it is simply the label of the nearest extended reference
  # node.  After grafting, the novel prototype sits at ~(8,0), so novel cells
  # must project to it.  A majority (>= 80%) must carry a "novel_" label.
  novel_direct_labels <- results$old_som_label[novel_idx]
  frac_novel_labelled <- mean(grepl("^novel_", novel_direct_labels))
  expect_gte(frac_novel_labelled, 0.8)
})


# ---------------------------------------------------------------------------
# Test 2: Per-batch artifact rejected; reproducible novel still minted
# ---------------------------------------------------------------------------
test_that("per-batch artifact is rejected while reproducible novel is minted", {
  withr::local_seed(2L)

  # Novel at (12, 0): very far from reference; score ~12 >> threshold.
  # Artifact at different coords per batch: (-8,5), (0,9), (8,5); score ~8.5-9.6.
  # In-dist at origin: score ~0.7.
  #
  # 3 groups x 200 in-dist + 100 novel + 100 artifact = 400/group = 1200 total.
  # tail_quantile=0.50: top 50% = 600 cells = all novel (300) + all artifact (300).
  # In-dist cells (~0.7) all fall below the threshold; threshold ~8.5.
  # Per group: 100 novel + 100 artifact = 200 tail cells.
  # min_cluster=70: k_g = min(8, 200, 200%/%70) = min(8, 200, 2) = 2.
  # Two clusters: novel at ~(12,0) with ~100 cells, artifact at ~(batch_loc) with ~100.
  # novel clusters link across groups (all near (12,0)) -> minted.
  # Artifact clusters are at (-8,5), (0,9), (8,5) -> pairwise distance > tol_factor*spacing.
  artifact_centers <- list(c(-8, 5), c(0, 9), c(8, 5))

  fix <- make_novelty_fixture(
    seed                 = 202L,
    n_groups             = 3L,
    n_ref_per_group      = 200L,
    novel_center         = c(12, 0),  # farther: score ~12, clearly > artifact ~8.5-9.6
    n_novel_per_group    = 100L,
    artifact_centers     = artifact_centers,
    n_artifact_per_group = 100L
  )

  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile = 0.50,   # top 50% = all novel + all artifact
    min_cluster   = 70L,
    min_batches   = 2L,
    tol_factor    = 2.0
  )

  expect_s3_class(cand, "somalign_novelty_candidates")

  # At least one candidate minted (the reproducible novel pop at (12,0))
  expect_gte(nrow(cand$prototypes), 1L)

  # The prototype closest to (12,0) must be within 2.0 units
  dists_to_novel <- apply(cand$prototypes, 1,
                           function(p) sqrt(sum((p - c(12, 0))^2)))
  expect_lt(min(dists_to_novel), 2.0)

  # Artifact centroids are in different groups only; confirm they did not mint.
  # Pairwise: (-8,5)-(0,9) = sqrt(64+16)=8.9; (-8,5)-(8,5) = 16; (0,9)-(8,5) = sqrt(64+16)=8.9.
  # All > tol_factor*node_spacing = 2*1 = 2.  None should be minted.
  for (ac in artifact_centers) {
    dists_to_ac <- apply(cand$prototypes, 1, function(p) sqrt(sum((p - ac)^2)))
    expect_gt(min(dists_to_ac), 3.0)
  }
})


# ---------------------------------------------------------------------------
# Test 3: min_batches respected (single-group novel NOT minted at min_batches=2)
# ---------------------------------------------------------------------------
test_that("population in only 1 group is not minted at min_batches = 2", {
  withr::local_seed(3L)

  # Base fixture: no cross-batch novel
  fix <- make_novelty_fixture(
    seed              = 303L,
    n_groups          = 3L,
    n_ref_per_group   = 200L,
    novel_center      = c(0, 0),  # won't be used
    n_novel_per_group = 0L
  )

  # Manually inject a single-batch novel population into batch1 only
  scaled_data <- fix$fit$query$scaled_data
  single_batch_cells <- matrix(
    c(rnorm(100, 10, 0.1),
      rnorm(100, 10, 0.1)),
    ncol = 2
  )
  colnames(single_batch_cells) <- colnames(scaled_data)

  # Augment query's scaled_data (patch the fit object)
  mock_fit <- fix$fit
  mock_fit$query$scaled_data <- rbind(scaled_data, single_batch_cells)
  groups_aug <- c(fix$groups, rep("batch1", 100))

  # top 15% of 700 cells = 105 cells; single-batch pop = 100 at (10,10).
  # They score high and enter the tail. k_g = min(8, 100, 100%/%70) = 1.
  # Cluster of 100 cells in batch1 only -> n_batches=1 < min_batches=2 -> NOT minted.
  cand <- somalign_novelty_candidates(
    mock_fit,
    groups_aug,
    tail_quantile = 0.85,
    min_cluster   = 70L,
    min_batches   = 2L,
    tol_factor    = 2.0
  )

  expect_s3_class(cand, "somalign_novelty_candidates")

  # Nothing near (10,10) should be minted
  if (nrow(cand$prototypes) > 0) {
    dists_to_single <- apply(cand$prototypes, 1,
                              function(p) sqrt(sum((p - c(10, 10))^2)))
    expect_gt(min(dists_to_single), 3.0)
  }
})


# ---------------------------------------------------------------------------
# Test 4: Score ranking — novel cells score higher than in-distribution cells
# ---------------------------------------------------------------------------
test_that("continuous score ranks novel cells above in-distribution cells", {
  withr::local_seed(4L)

  fix <- make_novelty_fixture(
    seed              = 404L,
    n_groups          = 2L,
    n_ref_per_group   = 200L,
    novel_center      = c(8, 0),
    n_novel_per_group = 100L
  )

  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile = 0.65,
    min_cluster   = 60L,
    min_batches   = 2L
  )

  n_cells <- nrow(fix$fit$query$scaled_data)
  expect_length(cand$score, n_cells)
  expect_true(all(is.finite(cand$score)))
  expect_true(is.numeric(cand$score))

  # Novel cells (near x=8) score much higher than in-dist cells (near origin)
  novel_idx   <- which(abs(fix$query_rows[, "x"] - 8) < 1.0)
  in_dist_idx <- which(abs(fix$query_rows[, "x"]) < 0.5)

  expect_gt(median(cand$score[novel_idx]),
            median(cand$score[in_dist_idx]))
  expect_gt(mean(cand$score[novel_idx]),
            mean(cand$score[in_dist_idx]))
})


# ---------------------------------------------------------------------------
# Test 5: Structural guards
# ---------------------------------------------------------------------------
test_that("structural guards reject bad inputs", {
  withr::local_seed(5L)

  fix <- make_novelty_fixture(
    seed              = 505L,
    n_groups          = 2L,
    n_ref_per_group   = 100L,
    n_novel_per_group = 0L
  )

  n_cells <- nrow(fix$fit$query$scaled_data)

  # Wrong group length
  expect_error(
    somalign_novelty_candidates(fix$fit, fix$groups[seq_len(n_cells - 1)]),
    regexp = "group.*length|length.*group",
    ignore.case = TRUE
  )

  # Wrong precomputed score length
  expect_error(
    somalign_novelty_candidates(fix$fit, fix$groups,
                                score = rnorm(n_cells - 1)),
    regexp = "score.*length|length.*score",
    ignore.case = TRUE
  )

  # Non-finite precomputed score
  bad_score    <- rnorm(n_cells)
  bad_score[1] <- NA
  expect_error(
    somalign_novelty_candidates(fix$fit, fix$groups, score = bad_score),
    regexp = "finite",
    ignore.case = TRUE
  )

  # Not a somalign_fit
  expect_error(
    somalign_novelty_candidates(list(), fix$groups),
    regexp = "somalign_fit",
    ignore.case = TRUE
  )
})


# ---------------------------------------------------------------------------
# Test 5b: Empty-tail edge case returns 0-row prototypes without error
# ---------------------------------------------------------------------------
test_that("very high tail_quantile returns 0 candidates without error", {
  withr::local_seed(6L)

  fix <- make_novelty_fixture(
    seed              = 606L,
    n_groups          = 2L,
    n_ref_per_group   = 100L,
    n_novel_per_group = 0L
  )

  # tail_quantile=0.999 -> only 0.1% of cells in tail; no group can accumulate
  # min_cluster cells.
  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile = 0.999,
    min_cluster   = 50L,
    min_batches   = 2L
  )

  expect_s3_class(cand, "somalign_novelty_candidates")
  expect_equal(nrow(cand$prototypes), 0L)
  expect_equal(ncol(cand$prototypes), ncol(fix$fit$query$scaled_data))
  expect_equal(colnames(cand$prototypes), colnames(fix$fit$query$scaled_data))
  expect_length(cand$n_groups_support, 0L)
  expect_length(cand$size,             0L)

  # Print method must work and display "candidates minted: 0"
  expect_output(print(cand), "candidates minted: 0")
})


# ---------------------------------------------------------------------------
# Test 5c: Print method works for non-empty candidates
# ---------------------------------------------------------------------------
test_that("print method works for non-empty somalign_novelty_candidates", {
  withr::local_seed(7L)

  fix <- make_novelty_fixture(
    seed              = 707L,
    n_groups          = 3L,
    n_ref_per_group   = 200L,
    novel_center      = c(8, 0),
    n_novel_per_group = 100L
  )

  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile = 0.65,
    min_cluster   = 60L,
    min_batches   = 2L
  )

  expect_output(print(cand), "somalign_novelty_candidates")
  expect_output(print(cand), "candidates minted")
})


# ---------------------------------------------------------------------------
# Deduplication tests
# ---------------------------------------------------------------------------
#
# Design:
# - node_spacing ≈ 1.0 for the 2x2 grid (median NN distance among the 4 nodes
#   at corners of a 1x1 square near origin; see existing tests).
# - tol_factor = 1.5  -> cross-batch matching window: 1.5 * 1 = 1.5
# - merge_tol_factor = 3.0 -> dedup window:            3.0 * 1 = 3.0
# - Two populations at (8,0) and (8,2.5): distance = 2.5,
#     WITHIN merge window (2.5 <= 3.0) but OUTSIDE matching window (2.5 > 1.5).
#   They produce two separate cross-group components post-matching, but should
#   merge into ONE prototype during dedup.
# - Both pops in all 3 groups -> union of group sets = {batch1,batch2,batch3}.
#   n_groups_support of merged candidate = 3 (union), not 6 (sum).
# - DISABLE test with same fixture (merge_tol_factor=NULL) must yield 2 candidates.
# - KEEP-FAR: two pops at (8,0) and (8,5): distance = 5.0 > 3.0 -> stay as TWO.

# Helper: make a fixture with TWO novel populations
make_two_novel_fixture <- function(
    seed    = 42L,
    center1 = c(8, 0),
    center2 = c(8, 2.5),
    n_groups          = 3L,
    n_ref_per_group   = 200L,
    n_novel_per_group = 100L
) {
  withr::local_seed(seed)

  ref_codebook <- matrix(
    c(-0.5, -0.5,
       0.5, -0.5,
      -0.5,  0.5,
       0.5,  0.5),
    nrow = 4, byrow = TRUE
  )
  colnames(ref_codebook) <- c("x", "y")

  ref <- somalign_reference_from_nodes(
    codebook    = ref_codebook,
    features    = c("x", "y"),
    center      = c(x = 0, y = 0),
    scale       = c(x = 1, y = 1),
    node_masses = rep(0.25, 4),
    label_prob  = matrix(
      c(1, 0,
        1, 0,
        0, 1,
        0, 1),
      nrow = 4, byrow = TRUE,
      dimnames = list(NULL, c("A", "B"))
    ),
    distance_quantiles = matrix(
      rep(c(0.5, 1.0, 1.5, 2.0), 4),
      nrow = 4, byrow = TRUE,
      dimnames = list(NULL, c("50%", "90%", "95%", "99%"))
    )
  )

  groups_vec <- character(0)
  query_rows <- matrix(numeric(0), ncol = 2)
  colnames(query_rows) <- c("x", "y")

  for (g in seq_len(n_groups)) {
    grp_name <- paste0("batch", g)

    # In-distribution cells near origin
    in_dist <- matrix(
      c(rnorm(n_ref_per_group, 0, 0.3),
        rnorm(n_ref_per_group, 0, 0.3)),
      ncol = 2
    )
    colnames(in_dist) <- c("x", "y")
    query_rows <- rbind(query_rows, in_dist)
    groups_vec <- c(groups_vec, rep(grp_name, n_ref_per_group))

    # Novel population 1
    nov1 <- matrix(
      c(rnorm(n_novel_per_group, center1[1], 0.1),
        rnorm(n_novel_per_group, center1[2], 0.1)),
      ncol = 2
    )
    colnames(nov1) <- c("x", "y")
    query_rows <- rbind(query_rows, nov1)
    groups_vec <- c(groups_vec, rep(grp_name, n_novel_per_group))

    # Novel population 2
    nov2 <- matrix(
      c(rnorm(n_novel_per_group, center2[1], 0.1),
        rnorm(n_novel_per_group, center2[2], 0.1)),
      ncol = 2
    )
    colnames(nov2) <- c("x", "y")
    query_rows <- rbind(query_rows, nov2)
    groups_vec <- c(groups_vec, rep(grp_name, n_novel_per_group))
  }

  # SOM: reference nodes + the two novel nodes
  query_codebook <- rbind(
    ref_codebook,
    matrix(c(center1[1], center1[2]), nrow = 1,
           dimnames = list(NULL, c("x", "y"))),
    matrix(c(center2[1], center2[2]), nrow = 1,
           dimnames = list(NULL, c("x", "y")))
  )
  query_som <- make_som(query_codebook)

  qry <- somalign_query(query_rows, ref, som_query = query_som)
  fit <- somalign_fit(qry, ref, solver = "internal")

  list(fit = fit, groups = groups_vec,
       center1 = center1, center2 = center2,
       n_groups = n_groups, query_rows = query_rows)
}


# ---------------------------------------------------------------------------
# Test MERGE-CLOSE: near-duplicate prototypes merge into one
# ---------------------------------------------------------------------------
test_that("MERGE-CLOSE: two near-duplicate novel pops merge into one candidate", {
  # Two pops at (8,0) and (8,2.5): Euclidean distance = 2.5
  # node_spacing ≈ 1.0 -> merge_tol_factor=3.0 -> threshold = 3.0
  # 2.5 <= 3.0 -> should merge.
  # tol_factor=1.5 -> matching window = 1.5 < 2.5 -> they form two separate
  # cross-group components before dedup (verified by DISABLE test below).
  # Both pops in {batch1, batch2, batch3} -> union n_groups_support = 3 (not 6).

  fix <- make_two_novel_fixture(
    seed    = 801L,
    center1 = c(8, 0),
    center2 = c(8, 2.5)
  )

  cand <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile    = 0.50,   # top 50% = all novel cells; 600 / 1200 total
    min_cluster      = 70L,    # 200 novel tail cells/group / 70 -> k_g=2
    min_batches      = 2L,
    tol_factor       = 1.5,
    merge_tol_factor = 3.0
  )

  expect_s3_class(cand, "somalign_novelty_candidates")

  # Exactly ONE prototype after dedup
  expect_equal(nrow(cand$prototypes), 1L)

  # n_groups_support = UNION = 3 (not 6 = sum of two 3-group components)
  expect_equal(cand$n_groups_support[1L], 3L)

  # Merged size = sum of both populations' cells
  expect_gt(cand$size[1L], 0L)

  # Merged prototype coordinates should be near the midpoint (~(8, 1.25))
  proto <- cand$prototypes[1, ]
  expect_lt(abs(proto["x"] - 8),   1.5)
  expect_lt(abs(proto["y"] - 1.25), 1.5)

  # merge_tol_factor recorded in params
  expect_equal(cand$params$merge_tol_factor, 3.0)
})


# ---------------------------------------------------------------------------
# Test DISABLE: merge_tol_factor=NULL reproduces pre-dedup two-prototype result
# ---------------------------------------------------------------------------
test_that("DISABLE: merge_tol_factor=NULL gives pre-dedup two separate candidates", {
  # Same fixture as MERGE-CLOSE but dedup disabled.
  # The two components at (8,0) and (8,2.5) (distance 2.5 > tol_factor*1=1.5)
  # should remain as TWO separate minted candidates.

  fix <- make_two_novel_fixture(
    seed    = 801L,   # same seed -> identical data
    center1 = c(8, 0),
    center2 = c(8, 2.5)
  )

  cand_no_dedup <- somalign_novelty_candidates(
    fix$fit,
    fix$groups,
    tail_quantile    = 0.50,
    min_cluster      = 70L,
    min_batches      = 2L,
    tol_factor       = 1.5,
    merge_tol_factor = NULL   # DISABLE
  )

  expect_s3_class(cand_no_dedup, "somalign_novelty_candidates")

  # TWO prototypes without dedup
  expect_equal(nrow(cand_no_dedup$prototypes), 2L)

  # merge_tol_factor recorded as NULL in params
  expect_null(cand_no_dedup$params$merge_tol_factor)
})


# ---------------------------------------------------------------------------
# Test KEEP-FAR: well-separated prototypes are not merged
# ---------------------------------------------------------------------------
test_that("KEEP-FAR: two distant novel pops stay as two separate candidates", {
  # Pops at (8,0) and (8,5): distance = 5.0
  # merge window = 3.0 * 1.0 = 3.0; 5.0 > 3.0 -> no merge.

  fix_far <- make_two_novel_fixture(
    seed    = 802L,
    center1 = c(8, 0),
    center2 = c(8, 5)
  )

  cand_far <- somalign_novelty_candidates(
    fix_far$fit,
    fix_far$groups,
    tail_quantile    = 0.50,
    min_cluster      = 70L,
    min_batches      = 2L,
    tol_factor       = 1.5,
    merge_tol_factor = 3.0
  )

  expect_s3_class(cand_far, "somalign_novelty_candidates")

  # TWO prototypes must remain
  expect_equal(nrow(cand_far$prototypes), 2L)

  # Each prototype should be near one of the two true novel centers
  dists_to_c1 <- apply(cand_far$prototypes, 1,
                        function(p) sqrt(sum((p - c(8, 0))^2)))
  dists_to_c2 <- apply(cand_far$prototypes, 1,
                        function(p) sqrt(sum((p - c(8, 5))^2)))
  expect_lt(min(dists_to_c1), 1.5)
  expect_lt(min(dists_to_c2), 1.5)
})
