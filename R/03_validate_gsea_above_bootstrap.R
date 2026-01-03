#' Bootstrap validation of enrichment for genes above a threshold
#'
#' This function evaluates whether genes classified as "above threshold" in a
#' GSEA ranking are enriched compared to random gene sets of the same size
#' sampled from the whole transcriptome. The statistic used is the sum of
#' ranking scores for the "above threshold" genes, compared to a bootstrap
#' null distribution obtained by repeated random sampling.
#'
#' @param result A list returned by \code{gsea_single_vector()}, containing
#'   at least \code{ranking} (named numeric vector) and \code{hits}.
#' @param threshold Numeric threshold applied to the ranking score to define
#'   "above" genes.
#' @param n_boot Number of bootstrap samples (random gene sets) to generate.
#'   Default is 5000.
#' @param title Character string specifying the title of the output plot.
#' @param base_size Base font size for the ggplot2 theme.
#'
#' @details
#' The function compares the observed statistic (sum of ranking scores for
#' genes above threshold within the gene set) to a null distribution obtained
#' by repeatedly sampling random gene sets of identical size from the entire
#' ranked transcriptome. The empirical p-value corresponds to the probability
#' that a random gene set achieves a statistic greater than or equal to the
#' observed value.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{plot}: A ggplot2 object showing the bootstrap null distribution
#'     and the observed statistic.
#'   \item \code{observed}: The observed statistic.
#'   \item \code{null_distribution}: Numeric vector of bootstrap statistics.
#'   \item \code{p_value}: Empirical p-value.
#'   \item \code{n_above}: Number of genes above threshold in the gene set.
#' }
#'
#' @examples
#' \dontrun{
#' res <- gsea_single_vector(expr, geneset)
#' v <- validate_gsea_above_bootstrap(
#'   result = res,
#'   threshold = 0,
#'   n_boot = 5000,
#'   title = "Bootstrap validation"
#' )
#' print(v$plot)
#' v$p_value
#' }
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate
#' @importFrom ggplot2 labs theme_minimal
#' @export
validate_gsea_above_bootstrap <- function(result,
                                          threshold = 0,
                                          n_boot = 5000,
                                          title = "Bootstrap validation of enrichment",
                                          base_size = 14) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  ranking <- result$ranking
  hits    <- result$hits

  if (is.null(names(ranking))) {
    stop("`result$ranking` must be a named numeric vector.")
  }

  # Scores of genes from the geneset
  geneset_scores <- ranking[hits]

  # Genes above threshold
  above_idx <- geneset_scores > threshold
  if (!any(above_idx)) {
    stop("No genes are above the threshold.")
  }

  above_scores <- geneset_scores[above_idx]
  n_above <- length(above_scores)

  # Observed statistic
  observed <- sum(above_scores)

  # Bootstrap sampling
  all_scores <- ranking
  N <- length(all_scores)

  set.seed(123)   # FIX: no stats::set.seed
  null_dist <- replicate(n_boot, {
    idx <- sample.int(N, size = n_above, replace = FALSE)
    sum(all_scores[idx])
  })

  # Empirical p-value
  p_val <- mean(null_dist >= observed)

  # Plot
  df_null <- data.frame(stat = null_dist)

  p <- ggplot2::ggplot(df_null, ggplot2::aes(x = .data$stat)) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "grey80",
      color = "grey50"
    ) +
    ggplot2::geom_vline(
      xintercept = observed,
      color = "darkorange1",
      size = 1.2
    ) +
    ggplot2::annotate(
      "text",
      x = observed,
      y = Inf,
      label = sprintf("Observed = %.3g\np = %.3g", observed, p_val),
      vjust = 1.5,
      hjust = 1.1,
      size = base_size * 0.35,
      color = "darkorange1"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "n_above = ", n_above,
        " | threshold = ", threshold,
        " | n_boot = ", n_boot
      ),
      x = "Bootstrap statistic (sum of scores)",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = base_size)

  list(
    plot = p,
    observed = observed,
    null_distribution = null_dist,
    p_value = p_val,
    n_above = n_above
  )
}

