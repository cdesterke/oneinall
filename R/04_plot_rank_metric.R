#' Plot the ranking metric with threshold-aware gene set highlighting
#'
#' This function visualizes the ranking metric used in GSEA and highlights
#' genes from the gene set according to whether their ranking score is above
#' or below a user-defined threshold. It provides a complementary view to the
#' running enrichment score by showing how strongly each gene contributes to
#' the ranking metric.
#'
#' @param result A list returned by \code{gsea_single_vector()}, containing:
#'   \itemize{
#'     \item \code{ranking}: Numeric vector of ranking scores (sorted decreasingly).
#'     \item \code{hits}: Logical vector indicating which genes belong to the gene set.
#'     \item \code{running_ES}: Numeric vector of the running enrichment score.
#'   }
#' @param threshold Numeric threshold applied to the ranking score to classify
#'   gene set members as above or below threshold.
#' @param title Character string specifying the plot title.
#' @param base_size Base font size for the ggplot2 theme.
#'
#' @details
#' The plot displays the ranking metric across the transcriptome, with vertical
#' bars marking the positions of genes belonging to the gene set. Bars are
#' colored according to whether the gene's ranking score exceeds the specified
#' threshold. No vertical line for the ES peak is drawn in this version.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # Simulated example
#' expr <- rnorm(1000)
#' names(expr) <- paste0("Gene", 1:1000)
#' geneset <- sample(names(expr), 50)
#'
#' # Compute GSEA result
#' res <- gsea_single_vector(expr, geneset)
#'
#' # Plot threshold-aware ranking metric
#' p <- plot_rank_metric_thresholded(
#'   result = res,
#'   threshold = 0,
#'   title = "Rank metric (threshold-aware)"
#' )
#'
#' print(p)
#' }
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_segment
#' @importFrom ggplot2 scale_x_continuous expansion scale_color_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @export
plot_rank_metric_thresholded <- function(result,
                                         threshold = 0,
                                         title = "Rank metric (threshold-aware)",
                                         base_size = 14) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  ranking <- result$ranking
  hits    <- result$hits

  df <- data.frame(
    position   = seq_along(ranking),
    rank_score = ranking,
    hits       = hits
  )

  # Classification based on threshold
  df$interpretation <- ifelse(
    df$hits & df$rank_score > threshold,
    "above_threshold",
    ifelse(df$hits, "below_threshold", "non_hit")
  )

  # Plot (no vertical ES line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$position, y = .data$rank_score)) +

    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +

    ggplot2::geom_line(color = "grey40", size = 0.8) +

    ggplot2::geom_segment(
      data = df[df$hits, ],
      ggplot2::aes(
        x = .data$position,
        xend = .data$position,
        y = 0,
        yend = .data$rank_score,
        color = .data$interpretation
      ),
      size = 0.7
    ) +

    ggplot2::scale_color_manual(
      values = c(
        "above_threshold" = "darkorange1",
        "below_threshold" = "steelblue",
        "non_hit" = "grey70"
      ),
      breaks = c("above_threshold", "below_threshold"),
      labels = c("Above threshold", "Below threshold"),
      name = "Gene status"
    ) +

    ggplot2::labs(
      title = title,
      x = "Rank in transcriptome",
      y = "Rank metric"
    ) +

    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size * 0.9),
      legend.text  = ggplot2::element_text(size = base_size * 0.8)
    )

  p
}
