#' Differential co-expression analysis with CS-CORE
#'
#' Test for differential co-expression between two groups using a permutation test.
#' CS-CORE is run on each group separately to obtain co-expression estimates.
#' Group labels are then permuted \code{n_permu} times at the biological sample level
#' to build a null distribution of differences in co-expression estimates.
#' Two-sided permutation p-values are returned for each gene pair.
#' @seealso \href{https://changsubiostats.github.io/CS-CORE/articles/differential_coexpression.html}{Differential co-expression tutorial}
#'
#' @param object A Seurat object containing single-cell RNA-seq data,
#'   subsetted to a single cell type.
#' @param group_label A character string giving the name of the \code{object@meta.data} column
#'   that contains the group labels (e.g. \code{"Status"}).
#' @param sample_label A character string giving the name of the \code{object@meta.data} column
#'   that identifies biological samples / donors (e.g. \code{"Donor"}).
#'   Permutations are performed at the sample level to preserve within-sample cell correlations.
#' @param genes A character vector of gene names (length \eqn{p}) for which co-expression
#'   will be estimated.
#' @param group_levels Optional. A length-2 character vector specifying the two group labels
#'   and their order, e.g. \code{c("Healthy", "COVID")}.
#'   The observed difference is computed as group 2 minus group 1.
#'   Defaults to \code{sort(unique(groups))}.
#' @param n_permu Integer; number of permutations. Defaults to \code{100}.
#' @param seed Optional integer seed passed to \code{set.seed()} before selecting
#'   permutation combinations. Defaults to \code{NULL}.
#' @param n_cores Integer; number of cores for parallelizing the permutation loop via
#'   \code{parallel::mclapply}. Defaults to \code{1L} (sequential).
#'   Values greater than 1 are not supported on Windows.
#' @param verbose Logical; if \code{FALSE} (default), suppresses all output from the
#'   underlying \code{CSCORE} calls and prints only a progress bar for the permutation
#'   loop. If \code{TRUE}, prints detailed \code{[INFO]} messages from each step.
#' @param ... Additional arguments passed to \code{\link{CSCORE}}.
#' @importFrom utils combn txtProgressBar setTxtProgressBar
#'
#' @return A list with four elements:
#' \describe{
#'   \item{est_group1}{A \eqn{p \times p} matrix of co-expression estimates for group 1.}
#'   \item{est_group2}{A \eqn{p \times p} matrix of co-expression estimates for group 2.}
#'   \item{obs_diff}{A \eqn{p \times p} matrix of observed differences in co-expression
#'     (\code{est_group2 - est_group1}).}
#'   \item{p_value}{A \eqn{p \times p} matrix of two-sided permutation p-values.
#'     Diagonal entries are \code{NA}.}
#' }
#' @export
#'
#' @examples
#' # See a full example at:
#' # https://changsubiostats.github.io/CS-CORE/articles/differential_coexpression.html
#'
diff_CSCORE <- function(object, group_label, sample_label, genes,
                        group_levels = NULL,
                        n_permu = 100, seed = NULL,
                        n_cores = 1L,
                        verbose = FALSE,
                        ...) {
  if (!group_label %in% colnames(object@meta.data)) {
    stop(sprintf("[ERROR] '%s' is not a column in object@meta.data.", group_label))
  }
  if (!sample_label %in% colnames(object@meta.data)) {
    stop(sprintf("[ERROR] '%s' is not a column in object@meta.data.", sample_label))
  }

  groups  <- object@meta.data[[group_label]]
  samples <- object@meta.data[[sample_label]]

  if (is.null(group_levels)) group_levels <- sort(unique(groups))
  if (length(group_levels) != 2) {
    stop("[ERROR] Exactly two group levels are required for differential co-expression analysis.")
  }

  # Sample-level group membership
  Status_id_tab <- table(samples, groups)
  group2_ids <- rownames(Status_id_tab)[Status_id_tab[, group_levels[2]] > 0]
  all_ids    <- rownames(Status_id_tab)

  # Helper: run CSCORE, suppressing its output when verbose = FALSE
  run_CSCORE <- function(obj) {
    if (verbose) {
      CSCORE(obj, genes = genes, ...)
    } else {
      result <- NULL
      invisible(utils::capture.output(suppressMessages(
        result <- CSCORE(obj, genes = genes, ...)
      )))
      result
    }
  }

  # Observed co-expression and difference
  if (verbose) cat(sprintf("[INFO] Running CS-CORE on %s (%d cells).\n",
                            group_levels[1], sum(groups == group_levels[1])))
  res1 <- run_CSCORE(object[, groups == group_levels[1]])

  if (verbose) cat(sprintf("[INFO] Running CS-CORE on %s (%d cells).\n",
                            group_levels[2], sum(groups == group_levels[2])))
  res2 <- run_CSCORE(object[, groups == group_levels[2]])
  obs_diff <- res2$est - res1$est

  # Permutation test
  n_group2  <- length(group2_ids)
  all_comb  <- combn(all_ids, n_group2)
  if (!is.null(seed)) set.seed(seed)
  random_combs <- sample(ncol(all_comb), n_permu, replace = FALSE)

  run_one_permu <- function(i) {
    perm2_ids <- all_comb[, random_combs[i]]
    perm1_ids <- all_ids[!all_ids %in% perm2_ids]
    null_res1 <- run_CSCORE(object[, samples %in% perm1_ids])
    null_res2 <- run_CSCORE(object[, samples %in% perm2_ids])
    null_res2$est - null_res1$est
  }

  if (verbose) {
    cat(sprintf("[INFO] Starting permutation test: %d permutations on %d core(s).\n",
                n_permu, n_cores))
  }

  # Run permutation loop with a progress bar.
  # When n_cores > 1, permutations are processed in chunks of n_cores so that
  # mclapply can be used while still updating the progress bar between chunks.
  pb        <- utils::txtProgressBar(min = 0, max = n_permu, style = 3)
  diff_list <- vector("list", n_permu)
  chunks    <- split(seq_len(n_permu),
                     ceiling(seq_len(n_permu) / n_cores))
  for (chunk in chunks) {
    diff_list[chunk] <- parallel::mclapply(chunk, run_one_permu, mc.cores = n_cores)
    utils::setTxtProgressBar(pb, max(chunk))
  }
  close(pb)

  if (verbose) cat("[INFO] Permutation test complete.\n")

  # Assemble 3-D null array and compute two-sided p-values
  p        <- length(genes)
  diff_mat <- array(NA_real_, dim = c(p, p, n_permu))
  for (i in seq_len(n_permu)) diff_mat[, , i] <- diff_list[[i]]

  diff_centered <- sweep(diff_mat, c(1, 2), obs_diff, FUN = "-")
  pval1   <- apply(diff_centered > 0, c(1, 2), mean)
  pval2   <- apply(diff_centered < 0, c(1, 2), mean)
  p_value <- pmin(pval1, pval2) * 2
  diag(p_value) <- NA

  dimnames(obs_diff) <- list(genes, genes)
  dimnames(p_value)  <- list(genes, genes)

  list(
    est_group1 = res1$est,
    est_group2 = res2$est,
    obs_diff   = obs_diff,
    p_value    = p_value
  )
}
