test_that("label_guided=TRUE concentrates transport on concordant node pairs", {
  ref <- tiny_reference()
  # tiny_reference has 3 nodes: left(A), middle(ambiguous), right(B)
  # Build a query with 2 nodes matching same A/B taxonomy
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # Assign label_prob with A/B matching tiny_reference taxonomy
  query_obj$label_prob <- rbind(
    node1 = c(A = 0.95, B = 0.05),
    node2 = c(A = 0.05, B = 0.95)
  )

  fit_guided <- somalign_fit(query_obj, ref, label_guided = TRUE, epsilon = 0.1)
  fit_free   <- somalign_fit(query_obj, ref, label_guided = FALSE, epsilon = 0.1)

  expect_s3_class(fit_guided, "somalign_fit")

  # With label guidance, query node 1 (A-dominant) should send more mass to
  # reference node 1 (A-dominant, left), and query node 2 (B-dominant) should
  # send more mass to reference node 3 (B-dominant, right).
  plan_guided <- fit_guided$transport_plan
  plan_free   <- fit_free$transport_plan

  # Node 1 (A): concordant pair is ref node 1 (col 1)
  expect_gt(plan_guided[1, 1], plan_free[1, 1] * 0.5)
  # Node 2 (B): concordant pair is ref node 3 (col 3)
  expect_gt(plan_guided[2, 3], plan_free[2, 3] * 0.5)

  # Discordant pairs receive penalty: node1->ref3 and node2->ref1 should shrink
  expect_lt(plan_guided[1, 3], plan_free[1, 3] + 1e-6)
  expect_lt(plan_guided[2, 1], plan_free[2, 1] + 1e-6)
})

test_that("label_guided=FALSE (default) gives same result as not passing the arg", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  query_obj$label_prob <- rbind(
    node1 = c(A = 0.95, B = 0.05),
    node2 = c(A = 0.05, B = 0.95)
  )

  fit_default <- somalign_fit(query_obj, ref, epsilon = 0.1)
  fit_false   <- somalign_fit(query_obj, ref, label_guided = FALSE, epsilon = 0.1)

  expect_equal(fit_default$transport_plan, fit_false$transport_plan)
})

test_that("label_guided=TRUE with NULL query$label_prob raises an error", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # make_som produces a plain SOM, so label_prob will be NULL
  expect_null(query_obj$label_prob)

  expect_error(
    somalign_fit(query_obj, ref, label_guided = TRUE),
    "query\\$label_prob is NULL"
  )
})

test_that("fully-disjoint label columns raise a descriptive error", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # Set label_prob with column names disjoint from reference (A, B)
  query_obj$label_prob <- rbind(
    node1 = c(X = 0.9, Y = 0.1),
    node2 = c(X = 0.1, Y = 0.9)
  )

  expect_error(
    somalign_fit(query_obj, ref, label_guided = TRUE),
    "no shared labels found"
  )
})

test_that("unlabeled nodes (max prob < 0.5) are never penalized", {
  # Test .somalign_build_label_mask directly with 3-class taxonomy
  # so we can produce nodes with max < 0.5 (e.g., equal spread across 3 classes)
  q_lp <- rbind(
    c(0.8, 0.1, 0.1),  # dominant class 1 (labeled)
    c(0.33, 0.34, 0.33) # near-equal spread -> max < 0.5 (unlabeled)
  )
  colnames(q_lp) <- c("A", "B", "C")

  r_lp <- rbind(
    c(0.1, 0.8, 0.1),  # dominant class 2 (labeled)
    c(0.33, 0.33, 0.34) # near-equal spread -> max < 0.5 (unlabeled)
  )
  colnames(r_lp) <- c("A", "B", "C")

  mask <- .somalign_build_label_mask(q_lp, r_lp)

  # Labeled query node 1 (dom A) vs labeled ref node 1 (dom B): discordant -> TRUE
  expect_true(mask[1, 1])
  # Unlabeled query node 2: entire row should be FALSE
  expect_true(all(mask[2, ] == FALSE))
  # Unlabeled ref node 2: entire col should be FALSE
  expect_true(all(mask[, 2] == FALSE))
})

test_that("somalign_fit_two_pass with label_guided=TRUE runs without error", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1.1, 0, 1.1, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, c("a", "b"))),
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  query_obj$label_prob <- rbind(
    node1 = c(A = 0.9, B = 0.1),
    node2 = c(A = 0.1, B = 0.9)
  )

  fit2 <- expect_no_error(
    somalign_fit_two_pass(query_obj, ref, label_guided = TRUE,
                          epsilon_global = 0.3, epsilon_local = 0.1)
  )
  expect_s3_class(fit2, "somalign_fit")
  expect_true(!is.null(fit2$two_pass))
})
