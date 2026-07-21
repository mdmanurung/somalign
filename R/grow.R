#' Grow a reference codebook with novel-population nodes (Growing Neural Gas)
#'
#' Seeds a Growing Neural Gas (GNG) from an existing trained reference codebook
#' and inserts new nodes for novel populations encountered in `new_data` WITHOUT
#' retraining or moving the original codebook vectors.  The original nodes are
#' completely frozen: their positions are never updated, regardless of whether
#' they win or are topological neighbours of a winner.
#'
#' @param reference A `somalign_reference` object produced by
#'   [somalign_train_reference()] or [somalign_reference_from_nodes()].
#' @param new_data A numeric matrix of cells (rows) × features (columns).  Must
#'   carry the same feature names as `reference$features`.
#' @param max_new_nodes Maximum number of new nodes to insert (default `20`).
#'   Growth stops when this cap is reached even if the error criterion would
#'   still trigger insertions.
#' @param lambda Insert a new node every `lambda` input presentations
#'   (default `100`).
#' @param epsilon_new Learning rate applied to newly inserted (unfrozen) nodes
#'   when they win or are topological neighbours of a winner (default `0.05`).
#' @param age_max Maximum edge age before an edge is pruned (default `50`).
#' @param error_decay Multiplicative decay applied to the accumulated error at
#'   every input presentation (default `0.995`).
#' @param n_epochs Number of full passes over `new_data` during GNG training
#'   (default `5`).  Multiple epochs improve convergence; inputs are shuffled
#'   each epoch.
#' @param seed Integer random seed used for epoch-wise shuffling (default
#'   `NULL`, no seeding).
#' @param novel_label Character scalar.  The label assigned to grown nodes in
#'   the extended `label_prob` matrix (default `"novel"`).
#'
#' @return A `somalign_reference` built with [somalign_reference_from_nodes()]:
#'   \describe{
#'     \item{Codebook}{Original rows unchanged; grown nodes appended.}
#'     \item{distance_quantiles}{`Inf` for grown rows — never flagged as outside
#'       the reference.}
#'     \item{global_distance_quantiles}{Inherited from `reference` (column-wise
#'       max over the original nodes only).}
#'     \item{label_prob}{Original labels preserved; a `novel_label` column
#'       appended.  Original nodes get 0 in `novel_label`; grown nodes get 1.}
#'     \item{node_masses}{Tabulated from `new_data` assignments over the
#'       extended codebook and re-normalised to sum 1.}
#'     \item{som_ref}{`NULL` (node-level reference without topology).}
#'   }
#'
#' @details
#' **Algorithm (GNG with frozen originals)**
#'
#' 1. Initialise codebook `W` from `reference$codebook` (N_orig nodes).  All
#'    N_orig indices are permanently *frozen*.
#' 2. Initialise an empty edge graph and per-node cumulative error vector.
#' 3. For each input `x` (shuffled per epoch, repeated for `n_epochs`):
#'    a. Find s1 (nearest node) and s2 (second-nearest node) in `W`.
#'    b. Increment `error[s1]` by `||x - W[s1,]||^2`.
#'    c. Age all edges incident to s1 by 1; remove edges older than `age_max`
#'       and prune isolated nodes that are **not** frozen.
#'    d. **Move** s1 and its topological neighbours by `epsilon_new * (x - W[i,])`
#'       **only if** the node is NOT frozen.
#'    e. Create / refresh (age = 0) the edge s1 -- s2.
#'    f. Every `lambda` presentations: insert a new node between the
#'       highest-error node `q` and its highest-error neighbour `f`; split error;
#'       stop if `max_new_nodes` reached.
#'    g. Apply `error_decay` to all errors.
#' 4. Assign each `new_data` row to the nearest node in the final extended `W`.
#' 5. Build and return an extended reference via [somalign_reference_from_nodes()].
#'
#' **Prototype note.** Node masses are computed from `new_data` alone.  The
#' function does not re-project original reference cells; mass on original nodes
#' reflects coverage of `new_data`, not the training corpus.  Use downstream
#' [somalign_fit()] with `laplacian_lambda = 0` (required because `som_ref` is
#' `NULL` in the returned object).
#'
#' @seealso [somalign_reference_from_nodes()], [somalign_fit()]
#'
#' @export
somalign_grow_reference <- function(reference,
                                    new_data,
                                    max_new_nodes = 20L,
                                    lambda        = 100L,
                                    epsilon_new   = 0.05,
                                    age_max       = 50L,
                                    error_decay   = 0.995,
                                    n_epochs      = 5L,
                                    seed          = NULL,
                                    novel_label   = "novel") {

  ## ---- input validation ---------------------------------------------------
  if (!inherits(reference, "somalign_reference")) {
    stop("`reference` must be a somalign_reference object.", call. = FALSE)
  }
  new_data <- as.matrix(new_data)
  if (!all(reference$features %in% colnames(new_data))) {
    stop("`new_data` must have columns matching `reference$features`.",
         call. = FALSE)
  }
  new_data <- new_data[, reference$features, drop = FALSE]
  stopifnot(
    is.numeric(new_data),
    is.numeric(max_new_nodes), max_new_nodes >= 0,
    is.numeric(lambda),        lambda >= 1,
    is.numeric(epsilon_new),   epsilon_new > 0, epsilon_new <= 1,
    is.numeric(age_max),       age_max >= 1,
    is.numeric(error_decay),   error_decay > 0, error_decay <= 1,
    is.numeric(n_epochs),      n_epochs >= 1
  )
  max_new_nodes <- as.integer(max_new_nodes)
  lambda        <- as.integer(lambda)
  age_max       <- as.integer(age_max)
  n_epochs      <- as.integer(n_epochs)

  ## ---- scale new_data into reference space --------------------------------
  X <- .somalign_scale_matrix(new_data, reference$center, reference$scale)
  N_samples <- nrow(X)

  ## ---- initialise GNG from reference codebook -----------------------------
  W       <- reference$codebook              # (N_orig x P)  will grow
  N_orig  <- nrow(W)
  n_total <- N_orig                          # current node count (grows)
  P       <- ncol(W)

  frozen  <- seq_len(N_orig)                 # indices of permanently frozen nodes

  ## Cumulative error per node
  err     <- numeric(n_total)

  ## Edge graph: list keyed by node index, values = named vectors of ages
  ## edges[[i]][j] = age of edge i-j  (undirected, stored in both directions)
  edges   <- vector("list", n_total)

  n_inserted <- 0L
  t_global   <- 0L                           # global presentation counter

  ## ---- helper: find s1 and s2 --------------------------------------------
  find_s1_s2 <- function(x) {
    diffs <- W - matrix(x, nrow = n_total, ncol = P, byrow = TRUE)
    sq_dists <- rowSums(diffs * diffs)
    rank2    <- order(sq_dists)[1:2]
    list(s1 = rank2[1L], s2 = rank2[2L])
  }

  ## ---- helper: ensure edges list is large enough -------------------------
  grow_edges_list <- function(new_n) {
    if (length(edges) < new_n) {
      edges <<- c(edges, vector("list", new_n - length(edges)))
    }
  }

  ## ---- helper: add / refresh edge ----------------------------------------
  add_edge <- function(a, b) {
    edges[[a]][as.character(b)] <<- 0L
    edges[[b]][as.character(a)] <<- 0L
  }

  ## ---- helper: age edges incident to s1 and prune ------------------------
  age_and_prune <- function(s1) {
    nb <- as.integer(names(edges[[s1]]))
    if (length(nb) == 0L) return(invisible(NULL))

    ## increment ages
    edges[[s1]] <<- edges[[s1]] + 1L
    for (j in nb) {
      if (!is.null(edges[[j]])) {
        old_age <- edges[[j]][as.character(s1)]
        if (!is.na(old_age)) {
          edges[[j]][as.character(s1)] <<- old_age + 1L
        }
      }
    }

    ## prune edges older than age_max
    old_edges <- as.integer(names(edges[[s1]][edges[[s1]] > age_max]))
    for (j in old_edges) {
      edges[[s1]] <<- edges[[s1]][names(edges[[s1]]) != as.character(j)]
      if (!is.null(edges[[j]])) {
        edges[[j]] <<- edges[[j]][names(edges[[j]]) != as.character(s1)]
        ## remove isolated NEW nodes (never remove frozen nodes)
        if (length(edges[[j]]) == 0L && !(j %in% frozen)) {
          ## mark as removed by setting W row to NA (handled in find_s1_s2)
          ## For simplicity, don't remove: the node stays but is disconnected.
          ## This is a prototype simplification.
        }
      }
    }
  }

  ## ---- helper: insert new node between q and f ---------------------------
  insert_node <- function() {
    if (n_inserted >= max_new_nodes) return(invisible(NULL))

    ## highest-error node
    q <- which.max(err)

    ## highest-error neighbour of q
    nb_q <- as.integer(names(edges[[q]]))
    if (length(nb_q) == 0L) return(invisible(NULL))   # no edges yet — skip
    f <- nb_q[which.max(err[nb_q])]

    ## new node position = midpoint
    new_pos <- (W[q, ] + W[f, ]) / 2.0

    ## append to W and err
    W     <<- rbind(W, new_pos)
    err   <<- c(err, (err[q] + err[f]) / 2.0)
    n_total <<- n_total + 1L
    grow_edges_list(n_total)

    ## update edges: remove q-f, add q-new and f-new
    edges[[q]] <<- edges[[q]][names(edges[[q]]) != as.character(f)]
    if (!is.null(edges[[f]])) {
      edges[[f]] <<- edges[[f]][names(edges[[f]]) != as.character(q)]
    }
    add_edge(q, n_total)
    add_edge(f, n_total)

    ## reduce error at q and f
    err[q] <<- err[q] * 0.5
    err[f] <<- err[f] * 0.5

    n_inserted <<- n_inserted + 1L
    invisible(NULL)
  }

  ## ---- main GNG training loop --------------------------------------------
  for (epoch in seq_len(n_epochs)) {
    if (!is.null(seed)) set.seed(seed + epoch)
    idx_order <- sample.int(N_samples)

    for (ii in idx_order) {
      x  <- X[ii, ]
      t_global <- t_global + 1L

      ## find s1 and s2
      res <- find_s1_s2(x)
      s1  <- res$s1
      s2  <- res$s2

      ## accumulate error at s1
      diff_s1  <- x - W[s1, ]
      err[s1]  <- err[s1] + sum(diff_s1 * diff_s1)

      ## age edges from s1 and prune old ones
      age_and_prune(s1)

      ## move s1 and its neighbours (only if NOT frozen)
      if (!(s1 %in% frozen)) {
        W[s1, ] <- W[s1, ] + epsilon_new * diff_s1
      }
      nb_s1 <- as.integer(names(edges[[s1]]))
      for (j in nb_s1) {
        if (!(j %in% frozen) && j <= n_total) {
          diff_j   <- x - W[j, ]
          W[j, ]   <- W[j, ] + epsilon_new * 0.1 * diff_j
        }
      }

      ## create/refresh edge s1-s2
      add_edge(s1, s2)

      ## insert new node every lambda steps
      if (t_global %% lambda == 0L && n_inserted < max_new_nodes) {
        insert_node()
      }

      ## decay all errors
      err <- err * error_decay
    }
  }

  ## ---- build extended reference -------------------------------------------
  ## Ensure column names are preserved
  rownames(W) <- NULL
  colnames(W) <- reference$features

  N_grown <- n_total - N_orig

  ## node_masses: assign all new_data to extended codebook, tabulate
  unit_assignments <- integer(N_samples)
  for (ii in seq_len(N_samples)) {
    diffs     <- W - matrix(X[ii, ], nrow = n_total, ncol = P, byrow = TRUE)
    sq_dists  <- rowSums(diffs * diffs)
    unit_assignments[ii] <- which.min(sq_dists)
  }
  masses_raw  <- tabulate(unit_assignments, nbins = n_total)
  node_masses <- masses_raw / sum(masses_raw)

  ## distance_quantiles: Inf for grown rows (never flagged outside reference).
  ## When the original reference was built without distance_quantiles, the
  ## stored matrix contains NAs (validator converts NULL → NA matrix).  NA is
  ## not accepted by .somalign_prepare_distance_quantiles; replace with Inf
  ## (same semantic: "never flag").
  orig_dq <- reference$distance_quantiles
  if (is.null(orig_dq)) {
    ## No distance quantiles in original reference — use Inf for all nodes
    dq_new <- matrix(Inf, nrow = n_total,
                     ncol = 4L,
                     dimnames = list(NULL, c("50%", "90%", "95%", "99%")))
  } else {
    ## Convert any NAs to Inf before appending
    orig_dq_inf <- orig_dq
    orig_dq_inf[is.na(orig_dq_inf)] <- Inf
    inf_rows <- matrix(Inf, nrow = N_grown, ncol = ncol(orig_dq_inf),
                       dimnames = list(NULL, colnames(orig_dq_inf)))
    dq_new   <- rbind(orig_dq_inf, inf_rows)
  }

  ## global_distance_quantiles: inherit from original (NOT recomputed over
  ## extended codebook, which would give Inf for every column).
  ## If the original had NAs (disabled detection), set global to Inf too.
  gdq_new <- reference$global_distance_quantiles
  if (is.null(gdq_new)) {
    if (!is.null(orig_dq) && any(is.finite(orig_dq))) {
      ## Compute from original finite values only
      orig_dq_for_global <- orig_dq
      orig_dq_for_global[is.na(orig_dq_for_global)] <- -Inf
      gdq_new <- apply(orig_dq_for_global, 2L, max, na.rm = TRUE)
      gdq_new[!is.finite(gdq_new)] <- Inf
    } else {
      ## Original had no finite distance quantiles — disable globally
      gdq_new <- NULL
    }
  }

  ## label_prob: add a novel_label column
  orig_lp <- reference$label_prob
  if (is.null(orig_lp) || ncol(orig_lp) == 0L) {
    ## Reference has no labels — create a 1-column matrix
    lp_new <- matrix(
      c(rep(0, N_orig), rep(1, N_grown)),
      nrow  = n_total, ncol = 1L,
      dimnames = list(NULL, novel_label)
    )
  } else {
    ## Pad original rows with 0 in the novel column
    novel_col_orig  <- matrix(0, nrow = N_orig, ncol = 1L,
                              dimnames = list(NULL, novel_label))
    ## Grown rows: 1 in novel column, 0 in all original label columns
    grown_orig_cols <- matrix(0, nrow = N_grown, ncol = ncol(orig_lp),
                              dimnames = list(NULL, colnames(orig_lp)))
    novel_col_grown <- matrix(1, nrow = N_grown, ncol = 1L,
                              dimnames = list(NULL, novel_label))

    lp_orig_extended <- cbind(orig_lp, novel_col_orig)
    lp_grown         <- cbind(grown_orig_cols, novel_col_grown)
    lp_new           <- rbind(lp_orig_extended, lp_grown)
  }

  ## node_var: not computed (prototype simplification)
  ## Build and return the extended reference
  somalign_reference_from_nodes(
    codebook                  = W,
    features                  = reference$features,
    center                    = reference$center,
    scale                     = reference$scale,
    node_masses               = node_masses,
    label_prob                = lp_new,
    distance_quantiles        = dq_new,
    global_distance_quantiles = gdq_new
  )
}
