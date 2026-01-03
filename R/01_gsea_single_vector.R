#' Single-sample GSEA (Gene Set Enrichment Analysis)
#'
#' Computes a single-sample GSEA enrichment score for a given expression vector
#' and a gene set. This function implements the classical running-sum statistic
#' described in Subramanian et al. (2005), including permutation-based
#' p-values and normalized enrichment score (NES).
#'
#' @param expr_col A named numeric vector of gene expression values.
#'   Gene names must be stored in \code{names(expr_col)}.
#' @param genes_input A character vector containing the genes of the gene set.
#' @param rank_method Ranking method used to transform the expression vector.
#'   One of \code{"zscore"}, \code{"log"}, or \code{"raw"}.
#' @param n_perm Number of permutations used to estimate the empirical p-value.
#'   Default is 5000.
#'
#' @details
#' The function ranks all genes in decreasing order and computes the running
#' enrichment score (ES) based on weighted hits and misses. Permutations are
#' performed by shuffling the hit positions while keeping the ranking fixed.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{ranking}: The ranked expression vector (decreasing order).
#'   \item \code{hits}: Logical vector indicating which ranked genes belong
#'     to the input gene set.
#'   \item \code{running_ES}: The running enrichment score curve.
#'   \item \code{ES}: The enrichment score.
#'   \item \code{NES}: The normalized enrichment score.
#'   \item \code{p_value}: Empirical p-value based on permutations.
#'   \item \code{genes_present}: Genes from the input gene set found in the
#'     expression vector.
#' }
#'
#' @examples
#' \dontrun{
#' expr <- rnorm(1000)
#' names(expr) <- paste0("Gene", 1:1000)
#' geneset <- sample(names(expr), 50)
#' res <- gsea_single_vector(expr, geneset)
#' }
#'
#' @export
gsea_single_vector <- function(expr_col,
                               genes_input,
                               rank_method = c("zscore", "log", "raw"),
                               n_perm = 5000) {

  # Validate ranking method
  rank_method <- match.arg(rank_method)

  # Genes present in the expression vector
  genes_present <- intersect(genes_input, names(expr_col))
  if (length(genes_present) == 0) {
    stop("No genes from the list are present in the vector of expression.")
  }

  # --- Ranking step ---
  if (rank_method == "zscore") {
    ranking <- as.numeric(scale(expr_col))   # FIX: no stats::scale
    names(ranking) <- names(expr_col)
  } else if (rank_method == "log") {
    ranking <- log2(expr_col + 1)
  } else {
    ranking <- expr_col
  }

  # Sort decreasing
  ranking_sorted <- sort(ranking, decreasing = TRUE)
  N <- length(ranking_sorted)

  # Identify hits
  hits <- names(ranking_sorted) %in% genes_present
  Nh <- sum(hits)
  Nm <- N - Nh

  if (Nh == 0) {
    stop("No genes from the gene set remain after sorting.")
  }

  # Standard GSEA weights (p = 1)
  rank_weights <- abs(ranking_sorted)

  # Normalization safeguard
  sum_hit_weights <- sum(rank_weights[hits])
  if (sum_hit_weights == 0) {
    stop("The weights of genes in the gene set are all zero.")
  }

  # Running ES
  Phit  <- cumsum(hits * rank_weights / sum_hit_weights)
  Pmiss <- cumsum((!hits) / Nm)
  running_ES <- Phit - Pmiss

  # Enrichment score
  ES <- running_ES[which.max(abs(running_ES))]

  # --- Permutation test ---
  set.seed(123)   # FIX: no stats::set.seed

  perm_ES <- replicate(n_perm, {
    perm_hits <- sample(hits, replace = FALSE)
    sum_perm_weights <- sum(rank_weights[perm_hits])

    if (sum_perm_weights == 0) {
      return(0)
    }

    Phit_p  <- cumsum(perm_hits * rank_weights / sum_perm_weights)
    Pmiss_p <- cumsum((!perm_hits) / Nm)
    running <- Phit_p - Pmiss_p
    running[which.max(abs(running))]
  })

  # Empirical p-value
  if (ES >= 0) {
    p_val <- mean(perm_ES >= ES)
  } else {
    p_val <- mean(perm_ES <= ES)
  }

  # Normalized ES
  NES <- ES / mean(abs(perm_ES))

  # Output
  list(
    ranking = ranking_sorted,
    hits = hits,
    running_ES = running_ES,
    ES = ES,
    NES = NES,
    p_value = p_val,
    genes_present = genes_present
  )
}

