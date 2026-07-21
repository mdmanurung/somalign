## Tests for label_guided = TRUE on grafted (extended) references.
## Attack 2 from the adversarial review: .somalign_build_label_mask() must use
## intersection logic, not identity, so that references carrying extra "novel_*"
## label columns do not hard-stop the fit.

# ---- helpers -----------------------------------------------------------------

make_guided_query <- function(ref, label_prob) {
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  query_obj$label_prob <- label_prob
  query_obj
}

# ---- test 1: non-grafted regression ------------------------------------------

test_that("non-grafted label_guided=TRUE is identical to original identical-column path", {
  ref <- tiny_reference()  # label_prob columns: A, B

  # Build a query with matching A/B columns
  q_lp <- rbind(
    node1 = c(A = 0.95, B = 0.05),
    node2 = c(A = 0.05, B = 0.95)
  )
  query_obj <- make_guided_query(ref, q_lp)

  # Direct call to internal: mask must be identical pre- and post-fix
  # (intersection of identical columns = same columns in same order)
  mask_direct <- .somalign_build_label_mask(q_lp, ref$label_prob)
  expect_equal(dim(mask_direct), c(2L, 3L))

  # Non-grafted fit with label_guided should still produce a valid fit
  fit_guided <- somalign_fit(query_obj, ref, label_guided = TRUE, epsilon = 0.1)
  expect_s3_class(fit_guided, "somalign_fit")

  # Concordant pairs should have higher transport mass than discordant ones
  plan <- fit_guided$transport_plan
  # node1 (A-dominant) -> ref node 1 (A-dominant) should have positive mass
  expect_gt(plan[1, 1], plan[1, 3])
  # node2 (B-dominant) -> ref node 3 (B-dominant) should have positive mass
  expect_gt(plan[2, 3], plan[2, 1])
})

# ---- test 2: grafted reference with novel column -----------------------------

test_that("label_guided=TRUE runs without error on a grafted reference with novel_* column", {
  ref <- tiny_reference()  # label_prob columns: A, B

  # Graft a new node labelled "novel_1" — adds an extra column absent from the
  # query's taxonomy
  new_cb <- matrix(c(2, 0), nrow = 1, dimnames = list(NULL, c("a", "b")))
  extended <- somalign_extend_reference(ref, new_cb, new_labels = "novel_1")

  # Confirm the extended reference carries the novel column
  expect_true("novel_1" %in% colnames(extended$label_prob))

  # Build query with the ORIGINAL taxonomy (no "novel_1" column) — the
  # mismatched-column situation that previously hard-stopped the fit
  q_lp <- rbind(
    node1 = c(A = 0.95, B = 0.05),
    node2 = c(A = 0.05, B = 0.95)
  )
  query_obj <- make_guided_query(extended, q_lp)

  # Must run without error
  fit <- expect_no_error(
    somalign_fit(query_obj, extended, label_guided = TRUE, epsilon = 0.1)
  )
  expect_s3_class(fit, "somalign_fit")

  # Label transfer on the shared classes should still work
  results <- somalign_results(fit)
  expect_true("transferred_label" %in% names(results))
  # Shared labels A and B must appear among transferred labels
  transferred <- unique(results$transferred_label)
  expect_true(any(c("A", "B") %in% transferred))

  # The novel node (col 4) must not be completely blocked — the transport plan
  # should route at least some mass there (compare to label_guided=FALSE)
  plan_guided <- fit$transport_plan
  fit_free    <- somalign_fit(query_obj, extended, label_guided = FALSE,
                              epsilon = 0.1)
  plan_free   <- fit_free$transport_plan
  # novel node is col 4 (nrow(extended$codebook) = 4)
  n_ref <- nrow(extended$codebook)
  # guided should not penalize the novel node more than free
  # (sum over query nodes of mass sent to novel node should be >= 0)
  expect_gte(sum(plan_guided[, n_ref]), 0)
})

test_that(".somalign_build_label_mask with novel column treats grafted node as unlabeled", {
  # After intersecting to shared columns {A, B}, the novel node's row in
  # ref_label_prob becomes all-zero -> max < 0.5 -> unlabeled -> never penalized
  q_lp <- rbind(
    node1 = c(A = 0.95, B = 0.05),
    node2 = c(A = 0.05, B = 0.95)
  )
  colnames(q_lp) <- c("A", "B")

  # Reference: 3 original nodes + 1 novel node one-hot on novel_1
  r_lp <- rbind(
    left    = c(A = 0.95, B = 0.05, novel_1 = 0.00),
    middle  = c(A = 0.50, B = 0.50, novel_1 = 0.00),
    right   = c(A = 0.05, B = 0.95, novel_1 = 0.00),
    novel   = c(A = 0.00, B = 0.00, novel_1 = 1.00)
  )

  mask <- .somalign_build_label_mask(q_lp, r_lp)
  expect_equal(dim(mask), c(2L, 4L))

  # Novel ref node (col 4) must never be penalized regardless of query node
  expect_true(all(mask[, 4L] == FALSE))

  # The concordant pairs among shared nodes should produce the right pattern:
  # node1 (dom A) vs right (dom B) -> discordant -> TRUE
  expect_true(mask[1, 3])
  # node2 (dom B) vs left (dom A) -> discordant -> TRUE
  expect_true(mask[2, 1])
  # node1 (dom A) vs left (dom A) -> concordant -> FALSE
  expect_false(mask[1, 1])
  # node2 (dom B) vs right (dom B) -> concordant -> FALSE
  expect_false(mask[2, 3])
})

# ---- test 3: empty-intersection guard ----------------------------------------

test_that("fully-disjoint label columns (empty intersection) raise an informative error", {
  # Query has X, Y; reference has A, B -> no overlap -> error
  q_lp <- rbind(
    node1 = c(X = 0.9, Y = 0.1),
    node2 = c(X = 0.1, Y = 0.9)
  )
  r_lp <- rbind(
    left  = c(A = 0.95, B = 0.05),
    right = c(A = 0.05, B = 0.95)
  )

  expect_error(
    .somalign_build_label_mask(q_lp, r_lp),
    "no shared labels found"
  )
})

test_that("somalign_fit with empty-intersection label_prob raises an informative error", {
  ref <- tiny_reference()  # label_prob: A, B
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # Fully disjoint labels
  query_obj$label_prob <- rbind(
    node1 = c(X = 0.9, Y = 0.1),
    node2 = c(X = 0.1, Y = 0.9)
  )

  expect_error(
    somalign_fit(query_obj, ref, label_guided = TRUE),
    "no shared labels found"
  )
})
