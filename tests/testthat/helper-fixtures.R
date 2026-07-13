make_anchored_fixture <- function(seed = 1L) {
  withr::local_seed(seed)
  p <- 3L
  ref_data <- rbind(
    matrix(rnorm(30 * p, mean = -3, sd = 0.5), ncol = p),
    matrix(rnorm(30 * p, mean =  3, sd = 0.5), ncol = p)
  )
  colnames(ref_data) <- paste0("F", seq_len(p))
  ref <- somalign_train_reference(
    ref_data, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  shift    <- rep(1.0, p)
  qry_data <- ref_data + matrix(shift, nrow = nrow(ref_data), ncol = p, byrow = TRUE)
  qry <- somalign_query(
    qry_data, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  anc_idx <- seq_len(20L)
  list(
    ref      = ref,
    qry      = qry,
    ref_data = ref_data,
    qry_data = qry_data,
    anc_idx  = anc_idx,
    anchor_old = ref_data[anc_idx, , drop = FALSE],
    anchor_new = qry_data[anc_idx, , drop = FALSE]
  )
}

make_som <- function(codebook) {
  codebook <- as.matrix(codebook)
  if (is.null(colnames(codebook)) && ncol(codebook) == 2) {
    colnames(codebook) <- c("a", "b")
  }
  structure(
    list(codes = list(data = codebook)),
    class = "kohonen"
  )
}

tiny_reference <- function() {
  codebook <- rbind(
    c(-1, 0),
    c(0, 0),
    c(1, 0)
  )
  colnames(codebook) <- c("a", "b")
  somalign_reference_from_nodes(
    codebook = codebook,
    features = colnames(codebook),
    center = c(a = 0, b = 0),
    scale = c(a = 1, b = 1),
    node_masses = c(0.4, 0.2, 0.4),
    label_prob = rbind(
      left = c(A = 0.95, B = 0.05),
      middle = c(A = 0.5, B = 0.5),
      right = c(A = 0.05, B = 0.95)
    ),
    distance_quantiles = matrix(
      c(0.2, 0.4, 0.6, 0.8,
        0.2, 0.4, 0.6, 0.8,
        0.2, 0.4, 0.6, 0.8),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("50%", "90%", "95%", "99%"))
    )
  )
}
