#' Create Pie Charts for TE Expression Classification
#'
#' Generates different types of charts showing the proportion of self-expressed 
#' vs. gene-dependent transposable elements across different TE classes (LINE,
#' SINE, LTR, DNA, etc.).
#'
#' @param res Data frame containing TE classification results. 
#' @param plot.title Character string. Title for the plot. If NULL (default),
#'   no title is displayed.
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "pdf", "svg", "eps", "png", "tiff", "jpeg". Default is "png".
#' @param width Numeric. Plot width in inches. Default is 28 (suitable for
#'   4-6 TE classes).
#' @param height Numeric. Plot height in inches. Default is 7.
#' @param output_folder Character string. Directory where plots will be saved.
#'   Default is current directory (".").
#' @param colors Named character vector. Colors for TE expression types.
#'   Default is c("dependent" = "#FF1F5B", "self" = "#009ADE"). Names must
#'   match values in \code{TE_expression} column.
#' @param labels Named character vector. Labels for TE expression types in
#'   legend. Default is c("dependent" = "Gene-dependent TEs",
#'   "self" = "Self-expressed TEs").
#'
#' @details
#' The function creates a pie chart for each TE class showing:
#' \itemize{
#'   \item **Self-expressed TEs**: TEs transcribed from their own promoter
#'   \item **Gene-dependent TEs**: TEs transcribed as part of gene runthrough
#' }
#'
#' Each pie slice is labeled with its percentage. TE classes are displayed
#' as separate facets, arranged horizontally (up to \code{max_cols} columns).
#'
#' @return Invisible NULL. Called for side effects (creating plot files).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'classified_TEs' has TE_expression and TE_class columns
#' TE_classify_pie(
#'   res = TE_results$res.TEs,
#'   plot.title = "TE Expression Classification",
#'   device = c("pdf", "png"),
#'   output_folder = "results/classification"
#' )
#' }
#'
#'   
TE_classify_pie <- function(res,
                            plot.title = NULL,
                            device = "png",
                            width = 28,
                            height = 7,
                            output_folder = ".",
                            colors = c("dependent" = "#FF1F5B", "self" = "#009ADE"),
                            labels = c("dependent" = "Gene-dependent TEs",
                                       "self" = "Self-expressed TEs"),
                            save = "all") {
  
  # Input validation
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res)) {
    stop("'res' must be a data frame.", call. = FALSE)
  }
  
  if (nrow(res) == 0L) {
    stop("'res' contains no data rows.", call. = FALSE)
  }
  
  # Check required columns
  required_cols <- c("TE_expression", "TE_class")
  missing_cols <- setdiff(required_cols, colnames(res))
  
  if (length(missing_cols) > 0L) {
    stop(
      "'res' must contain columns: ",
      paste(required_cols, collapse = ", "),
      "\nMissing: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!dir.exists(output_folder)) {
    tryCatch(
      dir.create(output_folder, recursive = TRUE),
      error = function(e) {
        stop(
          "Cannot create output directory '", output_folder, "': ",
          e$message,
          call. = FALSE
        )
      }
    )
  }
  
  
  # Create plot
  filename <- file.path(output_folder, paste0("pie_TE_classes_", save))
  
  p <- tryCatch(
    .create_pie_plot(res = res, 
                      plot.title = plot.title,
                      colors = colors,
                      labels = labels
                      ),
    error = function(e) {
      stop(
        "Failed to create pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot 
  tryCatch(
    printDevice(p, filename, device,
                width = width, height = height),
    error = function(e) {
      warning(
        "Failed to save pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  invisible(NULL)
}

#' Create Pie Plot for TE Classification
#'
#' Internal function to generate the actual pie chart visualization.
#'
#' @param res Data frame with TE_expression and TE_class columns
#' @param plot.title Character string for plot title
#' @param colors Named character vector of colors
#' @param labels Named character vector of legend labels
#'
#' @return A ggplot object
#' @keywords internal
#' @noRd
#'
.create_pie_plot <- function(res,
                             plot.title = NULL,
                             colors = c("dependent" = "#FF1F5B", "self" = "#009ADE"),
                             labels = c("dependent" = "Gene-dependent TEs",
                                        "self" = "Self-expressed TEs")) {

  # Calculate proportions by TE class
  res_summary <- res %>%
    dplyr::count(TE_expression, TE_class) %>%
    dplyr::group_by(TE_class) %>%
    dplyr::reframe(
      TE_class = TE_class,
      TE_expression = TE_expression,
      prop = n / sum(n) * 100
    )
  
  res_summary$TE_expression <- factor(res_summary$TE_expression)
  
  nTE_class <- length(unique(res$TE_class))
  
  p <- ggplot2::ggplot(
    res_summary,
    ggplot2::aes(x = "", y = .data$prop, fill = .data$TE_expression)
  ) + 
    ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme_void(base_size = 24) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$TE_class), 
      ncol = nTE_class
    )  +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%0.1f%%", prop)), 
                       position = ggplot2::position_stack(vjust=0.5),
                       size=8) +
    ggplot2::scale_fill_manual(
      labels = labels,
      values = colors
    ) + 
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 20)
    ) 

  # Add title if provided
  if (!is.null(plot.title)) {
    p <- p + ggplot2::ggtitle(plot.title)
  }
  
  p
}
