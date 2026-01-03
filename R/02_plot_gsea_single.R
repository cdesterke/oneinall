#' Plot a single-sample GSEA enrichment curve
#'
#' This function generates a classical GSEA-style enrichment plot from the output
#' of \code{gsea_single_vector()}. It displays the running enrichment score (ES),
#' the positions of gene hits, and highlights genes above or below a user-defined
#' threshold on the ranking metric.
#'
#' @param result A list returned by \code{gsea_single_vector()}, containing
#'   ranking, hits, running_ES, ES, NES, and p_value.
#' @param threshold Numeric value used to classify genes as above or below
#'   threshold based on their ranking score.
#' @param title Plot title. Default is \code{"GSEA-style enrichment plot"}.
#' @param base_size Base font size for the ggplot2 theme.
#'
#' @details
#' The function reproduces the classical GSEA visualization: a running enrichment
#' score curve, vertical bars marking the positions of genes in the gene set,
#' and a legend indicating how many genes fall above or below the threshold.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{plot}: A ggplot2 object containing the enrichment plot.
#'   \item \code{geneset_table}: A data frame containing only the genes from the
#'     input gene set, including their ranking position, ranking score, threshold
#'     classification (above/below), and bar height used in the plot.
#' }
#'
#' @examples
#' \dontrun{
#' res <- gsea_single_vector(expr, geneset)
#' p <- plot_gsea_single(res, threshold = 0)
#' print(p$plot)
#' head(p$geneset_table)
#' }
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_segment
#' @importFrom ggplot2 scale_x_continuous expansion scale_color_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @export
plot_gsea_single <- function(
    result,
    threshold = 0,
    title = "GSEA-style enrichment plot",
    base_size = 14
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  # Build dataframe
  df <- data.frame(
    gene         = names(result$ranking),
    position     = seq_along(result$running_ES),
    rank_score   = result$ranking,
    running_ES   = result$running_ES,
    hits         = result$hits
  )

  # Subset hits
  df_hits <- df[df$hits, ]

  # Classification above/below threshold
  df_hits$interpretation <- ifelse(
    df_hits$rank_score > threshold,
    "above_threshold",
    "below_threshold"
  )

  # Bar height for visualization
  df_hits$bar_height <- ifelse(
    df_hits$interpretation == "above_threshold",
    0.25, 0.12
  )

  # Counts for legend
  n_above <- sum(df_hits$interpretation == "above_threshold")
  n_below <- sum(df_hits$interpretation == "below_threshold")

  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$position, y = .data$running_ES)) +

    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +

    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +

    ggplot2::geom_segment(
      data = df_hits,
      ggplot2::aes(
        x = .data$position,
        xend = .data$position,
        y = 0,
        yend = .data$bar_height,
        color = .data$interpretation
      ),
      inherit.aes = FALSE,
      size = 0.7
    ) +

    ggplot2::scale_color_manual(
      name = "Gene status",
      values = c(
        "above_threshold" = "darkorange1",
        "below_threshold" = "seagreen3"
      ),
      labels = c(
        paste0("Above threshold (", n_above, ")"),
        paste0("Below threshold (", n_below, ")")
      )
    ) +

    ggplot2::geom_line(color = "steelblue", size = 1.2) +

    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "ES = ", round(result$ES, 3),
        " | NES = ", round(result$NES, 3),
        " | p = ", signif(result$p_value, 3),
        " | threshold (rank) = ", signif(threshold, 3)
      ),
      x = "Rank in transcriptome",
      y = "Running Enrichment Score"
    ) +

    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size * 0.9),
      legend.text  = ggplot2::element_text(size = base_size * 0.8)
    )

  list(
    plot = p,
    geneset_table = df_hits
  )
}
