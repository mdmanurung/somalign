## Anchor-benefit validation (does using repeat samples improve label transfer).

skip_if_not_installed("kohonen")

# A labelled two-batch fixture with repeat (anchor) samples and a batch shift.
make_anchor_benefit_fixture <- function(seed = 1L, shift = 2.0) {
  withr::local_seed(seed)
  p <- 3L
  n_per <- 60L
  # reference batch: two well-separated populations, labelled
  ref_data <- rbind(
    matrix(rnorm(n_per * p, -3, 0.5), ncol = p),
    matrix(rnorm(n_per * p,  3, 0.5), ncol = p)
  )
  colnames(ref_data) <- paste0("F", seq_len(p))
  ref_labels <- rep(c("A", "B"), each = n_per)
  ref <- somalign_train_reference(ref_data, labels = ref_labels,
                                  grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 20)
  # query batch: same biology, shifted by `shift` in every marker
  qshift <- matrix(shift, nrow = nrow(ref_data), ncol = p, byrow = TRUE)
  qry_data <- ref_data + qshift
  qry_labels <- ref_labels
  qry <- somalign_query(qry_data, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 20)
  # anchors: repeat cells, unshifted in ref space and shifted in query space
  anc_idx <- c(1:15, (n_per + 1):(n_per + 15))
  list(ref = ref, qry = qry, qry_labels = qry_labels,
       anchor_old = ref_data[anc_idx, , drop = FALSE],
       anchor_new = ref_data[anc_idx, , drop = FALSE] +
         matrix(shift, nrow = length(anc_idx), ncol = p, byrow = TRUE))
}

test_that("anchor_benefit grid has one row per rho with sane metrics", {
  fx <- make_anchor_benefit_fixture()
  ab <- somalign_anchor_benefit(
    fx$qry, fx$ref, fx$qry_labels, fx$anchor_old, fx$anchor_new,
    rho_grid = c(0, 1, 5), epsilon = 0.1
  )
  expect_s3_class(ab, "somalign_anchor_benefit")
  expect_equal(nrow(ab$grid), 3L)
  expect_true(all(c("rho_anchor", "accuracy", "macro_f1", "mcc", "coverage", "ece")
                  %in% names(ab$grid)))
  expect_true(all(ab$grid$accuracy >= 0 & ab$grid$accuracy <= 1))
  expect_true(all(ab$grid$coverage >= 0 & ab$grid$coverage <= 1))
  expect_equal(nrow(ab$baseline), 1L)
  expect_equal(ab$baseline$rho_anchor, 0)
  op <- withVisible(print(ab))
  expect_false(op$visible)
})

test_that("rho_anchor = 0 matches plain somalign_fit label metrics", {
  fx <- make_anchor_benefit_fixture()
  ab <- somalign_anchor_benefit(
    fx$qry, fx$ref, fx$qry_labels, fx$anchor_old, fx$anchor_new,
    rho_grid = 0, epsilon = 0.1, solver = "internal"
  )
  fit <- somalign_fit(fx$qry, fx$ref, epsilon = 0.1, solver = "internal")
  res <- somalign_results(fit, include_correction = FALSE)
  m <- somalign_label_metrics(res$transferred_label, fx$qry_labels,
                              res$transferred_label_accepted)
  expect_equal(ab$grid$accuracy[1], m$accuracy, tolerance = 1e-8)
  expect_equal(ab$grid$mcc[1], m$mcc, tolerance = 1e-8)
  expect_equal(ab$grid$coverage[1], m$coverage, tolerance = 1e-8)
})

test_that("efficient path equals full somalign_fit_anchored at a positive rho", {
  fx <- make_anchor_benefit_fixture()
  rho <- 2
  ab <- somalign_anchor_benefit(
    fx$qry, fx$ref, fx$qry_labels, fx$anchor_old, fx$anchor_new,
    rho_grid = rho, epsilon = 0.1, solver = "internal"
  )
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
    rho_anchor = rho, correction = "cost_bonus", epsilon = 0.1, solver = "internal"
  )
  res <- somalign_results(fit, include_correction = FALSE)
  m <- somalign_label_metrics(res$transferred_label, fx$qry_labels,
                              res$transferred_label_accepted)
  # Same plan -> same labels -> identical metrics.
  expect_equal(ab$grid$accuracy[1], m$accuracy, tolerance = 1e-8)
  expect_equal(ab$grid$macro_f1[1], m$macro_f1, tolerance = 1e-8)
  expect_equal(ab$grid$mcc[1], m$mcc, tolerance = 1e-8)
})

test_that("eval_mask restricts scoring to selected cells", {
  fx <- make_anchor_benefit_fixture()
  mask <- rep(c(TRUE, FALSE), length.out = nrow(fx$qry$scaled_data))
  ab <- somalign_anchor_benefit(
    fx$qry, fx$ref, fx$qry_labels, fx$anchor_old, fx$anchor_new,
    rho_grid = 0, eval_mask = mask
  )
  # scored count equals accepted cells within the mask
  fit <- somalign_fit(fx$qry, fx$ref, epsilon = 0.1, solver = "internal")
  res <- somalign_results(fit, include_correction = FALSE)
  m <- somalign_label_metrics(res$transferred_label[mask], fx$qry_labels[mask],
                              res$transferred_label_accepted[mask])
  expect_equal(ab$grid$accuracy[1], m$accuracy, tolerance = 1e-8)
})

test_that("anchor_benefit validates inputs", {
  fx <- make_anchor_benefit_fixture()
  expect_error(
    somalign_anchor_benefit(fx$qry, fx$ref, fx$qry_labels[1:5],
                            fx$anchor_old, fx$anchor_new),
    "one entry per query cell"
  )
  expect_error(
    somalign_anchor_benefit(fx$qry, fx$ref, fx$qry_labels,
                            fx$anchor_old, fx$anchor_new, rho_grid = c(-1, 1)),
    "non-negative"
  )
})
