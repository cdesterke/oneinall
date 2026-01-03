#' GSEA-style enrichment plot with leading-edge highlighting
#'
#' This function generates a full GSEA-style enrichment plot including the
#' running enrichment score (ES), hit positions, and automatic extraction and
#' highlighting of the leading-edge subset (genes contributing to the maximum
#' absolute ES).
#'
#' @param result A list returned by \code{gsea_single_vector()}, containing
#'   \code{ranking}, \code{hits}, \code{running_ES}, \code{ES}, \code{NES},
#'   and \code{p_value}.
#' @param title Character string specifying the plot title.
#' @param base_size Base font size for the ggplot2 theme.
#'
#' @details
#' The leading-edge subset corresponds to the genes in the gene set that appear
#' before the position of the maximum absolute running enrichment score. These
#' genes are highlighted on the enrichment plot.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{plot}: A ggplot2 object showing the GSEA-style plot.
#'   \item \code{leading_edge}: A data frame of leading-edge genes.
#' }
#'
#' @examples
#' \dontrun{
#' res <- gsea_single_vector(expr, geneset)
#' p <- plot_gsea_with_leading_edge(res, title = "random genes in cells", base_size = 16)
#' print(p$plot)
#' head(p$leading_edge)
#' }
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_segment geom_vline
#' @importFrom ggplot2 scale_x_continuous expansion labs theme_minimal theme element_text
#' @export
plot_gsea_with_leading_edge <- function(result,
                                        title = "GSEA-style enrichment plot",
                                        base_size = 14) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  ranking <- result$ranking
  hits    <- result$hits
  running <- result$running_ES

  df <- data.frame(
    gene       = names(ranking),
    position   = seq_along(ranking),
    rank_score = ranking,
    running_ES = running,
    hits       = hits,
    stringsAsFactors = FALSE
  )

  # --- Leading edge extraction ---
  max_pos <- which.max(abs(running))
  leading_edge <- df[df$hits & df$position <= max_pos, ]

  # --- Plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$position, y = .data$running_ES)) +

    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.01))) +

    ggplot2::geom_line(color = "steelblue", size = 1.2) +

    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +

    ggplot2::geom_segment(
      data = df[df$hits, ],
      ggplot2::aes(
        x = .data$position,
        xend = .data$position,
        y = 0,
        yend = 0.15
      ),
      color = "grey40",
      size = 0.6
    ) +

    ggplot2::geom_segment(
      data = leading_edge,
      ggplot2::aes(
        x = .data$position,
        xend = .data$position,
        y = 0,
        yend = 0.25
      ),
      color = "darkorange1",
      size = 0.9
    ) +

    ggplot2::geom_vline(
      xintercept = max_pos,
      color = "red",
      linetype = "dashed",
      size = 0.9
    ) +

    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "ES = ", round(result$ES, 3),
        " | NES = ", round(result$NES, 3),
        " | p = ", signif(result$p_value, 3),
        " | Leading edge = ", nrow(leading_edge), " genes"
      ),
      x = "Rank in transcriptome",
      y = "Running Enrichment Score"
    ) +

    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  list(
    plot = p,
    leading_edge = leading_edge
  )
}
