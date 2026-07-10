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
