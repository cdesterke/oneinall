#' Extract leading-edge genes from a GSEA result
#'
#' This function extracts the leading-edge subset from the output of
#' \code{gsea_single_vector()}. The leading-edge genes are those that contribute
#' to the maximum absolute enrichment score (ES), following the definition used
#' in the original GSEA algorithm.
#'
#' @param result A list returned by \code{gsea_single_vector()}, containing
#'   \code{ranking}, \code{hits}, and \code{running_ES}.
#'
#' @details
#' The leading-edge subset corresponds to the genes in the gene set that appear
#' before the position of the maximum absolute running enrichment score. These
#' genes are considered the main contributors to the enrichment signal.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{gene}: Gene names.
#'   \item \code{position}: Rank position in the transcriptome.
#'   \item \code{rank_score}: Ranking metric value.
#'   \item \code{running_ES}: Running enrichment score at that position.
#' }
#'
#' @examples
#' \dontrun{
#' res <- gsea_single_vector(expr, geneset)
#' le <- extract_leading_edge(res)
#' head(le)
#' }
#'
#' @export
extract_leading_edge <- function(result) {

  # Basic structure checks
  if (is.null(result$ranking) || is.null(result$hits) || is.null(result$running_ES)) {
    stop("`result` must contain ranking, hits, and running_ES.")
  }

  ranking <- result$ranking
  hits    <- result$hits
  running <- result$running_ES

  if (length(ranking) != length(hits) ||
      length(ranking) != length(running)) {
    stop("`ranking`, `hits`, and `running_ES` must have the same length.")
  }

  # Position of maximum absolute ES
  max_pos <- which.max(abs(running))

  # Build full data frame
  df <- data.frame(
    gene       = names(ranking),
    position   = seq_along(ranking),
    rank_score = ranking,
    running_ES = running,
    hits       = hits,
    stringsAsFactors = FALSE
  )

  # Leading-edge = hits before the ES peak
  leading_edge <- df[df$hits & df$position <= max_pos, ]

  # Remove hits column for clarity
  leading_edge$hits <- NULL

  leading_edge
}
