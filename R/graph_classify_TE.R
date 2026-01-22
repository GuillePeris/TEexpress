#' Create Classification Plots for TE Expression
#'
#' Generates multiple types of plots showing the proportion of self-expressed
#' vs. gene-dependent transposable elements across different TE classes (LINE,
#' SINE, LTR, DNA, etc.).
#'
#' @param res Data frame containing TE classification results. Must include
#'   columns: \code{TE_expression} (classification: "dependent" or "self"),
#'   and \code{expression_type} (detailed classification for stacked bar plots).
#' @param plot.title Character string. Title for the plot. If NULL (default),
#'   no title is displayed.
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "pdf", "svg", "eps", "png", "tiff", "jpeg". Default is "png".
#' @param width.pie Numeric. Plot width in inches for pie graph. Default is 14 (suitable for
#'   4-6 TE classes).
#' @param height.pie Numeric. Plot height in inches for pie graph. Default is 7.
#' @param height.bar Numeric. Plot height in inches for stack bar. Default is 7.
#' @param width.bar Numeric. Plot width in inches for stack. Default is 7 (suitable for
#'   4-6 TE classes).
#' @param output_folder Character string. Directory where plots will be saved.
#'   Default is current directory (".").
#' @param colors Named character vector. Colors for TE expression types.
#'   Default is c("dependent" = "#FF1F5B", "self" = "#009ADE"). Names must
#'   match values in \code{TE_expression} column.
#' @param labels Named character vector. Labels for TE expression types in
#'   legend. Default is c("dependent" = "Gene-dependent TEs",
#'   "self" = "Self-expressed TEs").
#' @param save Character string. Which TEs to save in output file:
#'   \itemize{
#'     \item "all": All expressed TEs (default)
#'     \item "dys": Only significantly dysregulated TEs (up or down)
#'     \item "up": Only significantly upregulated TEs
#'     \item "down": Only significantly downregulated TEs
#'   }
#'
#' The function creates two types of plots:
#'
#' \strong{1. Pie Charts:}
#' One pie chart per TE class showing the overall proportion of:
#' \itemize{
#'   \item \strong{Self-expressed TEs}: Transcribed from their own promoter
#'   \item \strong{Gene-dependent TEs}: Transcribed as part of gene runthrough
#' }
#'
#' \strong{2. Stacked Bar Charts:}
#' Grouped bars showing detailed classification breakdown:
#' \itemize{
#'   \item Exon, Intron (expressed), Intron (not expressed), etc.
#'   \item Grouped by self-expressed vs. gene-dependent
#'   \item One panel per TE class
#' }
#'
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
                            width.pie = 14,
                            height.pie = 7,
                            width.bar = 7,
                            height.bar = 7,
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
  
  # Check required elements
  required_columns <- c("expression_type", "TE_expression", 
                        "TE_name", "TE_family", "TE_class")
  missing_columns <- setdiff(required_columns, colnames(res))
  
  if (length(missing_columns) > 0L) {
    stop(
      "'res' must contain columns: ",
      paste(required_columns, collapse = ", "),
      "\nMissing: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ============================================================
  # Create and Save Pie Charts
  # ============================================================  
  
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
                width = width.pie, height = height.pie),
    error = function(e) {
      warning(
        "Failed to save pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # ============================================================
  # Create and Save Stacked Bar Charts
  # ============================================================
  
  filename <- file.path(output_folder, paste0("stack_bar_TE_classes_", save))
  
  p <- tryCatch(
    .create_grouped_stack_bar(res = res, 
                              plot.title = plot.title
                             ),
    error = function(e) {
      stop(
        "Failed to create stack bar TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot 
  tryCatch(
    printDevice(p, filename, device,
                width = width.bar, height = height.bar),
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
    dplyr::count(.data$TE_expression, .data$TE_class) %>%
    dplyr::group_by(.data$TE_class) %>%
    dplyr::reframe(
      TE_class = .data$TE_class,
      TE_expression = .data$TE_expression,
      prop = .data$n / sum(.data$n) * 100
    )
  
  if (nrow(res_summary) == 0L) {
    stop("No data remaining after summarization.", call. = FALSE)
  }
  
  # Ensure factor for consistent ordering
  res_summary$TE_expression <- factor(
    res_summary$TE_expression,
    levels = c("dependent", "self")
  )
  
  nTE_class <- length(unique(res_summary$TE_class))
  
  # Create pie charts
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
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%0.1f%%", .data$prop)), 
                       position = ggplot2::position_stack(vjust=0.5),
                       size=6) +
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

#' Create Grouped Stacked Bar Chart for TE Classification
#'
#' Internal function to generate stacked bar charts showing detailed
#' classification breakdown (Exon, Intron expressed/not expressed, etc.)
#' grouped by self-expressed vs. gene-dependent.
#'
#' @param res Data frame with TE_expression, expression_type, and TE_class columns
#' @param plot.title Character string for plot title
#'
#' @return A ggplot object
#' @keywords internal
#' @noRd
.create_grouped_stack_bar <- function(res,
                                      plot.title = NULL) {
 
  # Calculate proportions by TE class
  res_summary <- res %>%
    dplyr::count(.data$TE_expression, 
                 .data$expression_type, 
                 .data$TE_class
    ) %>%
    dplyr::group_by(.data$TE_class) %>%
    dplyr::reframe(
      TE_class = .data$TE_class,
      TE_expression = .data$TE_expression,
      expression_type = .data$expression_type,
      prop = .data$n / sum(.data$n) * 100
    )
  
  if (nrow(res_summary) == 0L) {
    stop("No data remaining after summarization.", call. = FALSE)
  }
  
  # Ensure factor for consistent ordering
  res_summary$TE_expression <- factor(
    res_summary$TE_expression,
    levels = c("dependent", "self")
  )
  res_summary$expression_type <- factor(res_summary$expression_type)
  
  # Shorten label for better display
  levels(res_summary$TE_expression)[
    levels(res_summary$TE_expression) == "dependent"
  ] <- "dep."  
  
  # Create stacked bar chart
  p <- ggplot2::ggplot(
    res_summary,
    ggplot2::aes(x = .data$TE_expression, 
                 y = .data$prop, 
                 fill = .data$expression_type,
                 ymin = 0,
                 ymax = max(.data$prop + 5, 100)
    )
  ) + 
    ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$TE_class)
  )   
  
  
  # Add percentage labels
  p <- p +
    ggfittext::geom_fit_text(
      ggplot2::aes(label = sprintf("%0.1f%%", .data$prop)),
      min.size = 8,
      size = 16,
      position = ggplot2::position_stack(vjust = 0.5),
      contrast = TRUE,
      show.legend = FALSE
    )
  
  # Apply theme
  p <- p +
    ggplot2::theme(
      # Grid
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),

      # Background
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Axis
      axis.text = ggplot2::element_text(colour = "grey30", size = 12),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        size = 20,
        colour = "grey30",
        hjust = 0.5,
        margin = ggplot2::margin(r = 10)
      ),
      axis.ticks.y = ggplot2::element_line(color = "grey30", linewidth = 0.5),
      axis.ticks.length.y = ggplot2::unit(0.25, "cm"),
      
      # Legend
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 24, color = "red")
    ) +
    ggplot2::labs(
      y = "% TE loci"
    )
  
  # Format y-axis as percentages
  p <- p +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, " %"),
      expand = c(0, 1)
    )
  
  
  # Add title if provided
  if (!is.null(plot.title)) {
    p <- p + ggplot2::ggtitle(plot.title)
  }
  
  p 
}

#' Create Donut Pie Chart for Single TE Type
#'
#' Creates a donut/pie chart for a specific TE type showing the breakdown
#' of expression types. Requires the \code{webr} package.
#' 
#' Please, note that \code{webr} package is not maintained and this function
#' throws several warnings of deprecated uses.
#'
#' @param res Data frame with classification results
#' @param TE_feature Character string. Specific TE to plot (e.g., "LINE", "L1")
#' @param plot.title Character string. Plot title
#' @param output_folder Character string. Output directory
#' @param width Numeric. Plot width in inches. Default is 7.
#' @param height Numeric. Plot height in inches. Default is 7.
#' @param device Character vector. Output format(s). Default is "png".
#' @param prefix Character string. Optional prefix for filename
#'
#' @return Invisible NULL. Called for side effects (creating plot file).
#'
#' @export
#' #' @examples
#' \dontrun{
#' # Requires webr package
#' if (requireNamespace("webr", quietly = TRUE)) {
#'   create_pie_donut(
#'     res = classified_TEs,
#'     TE_feature = "LINE",
#'     output_folder = "plots",
#'     prefix = "detailed"
#'   )
#' }
#' }
create_pie_donut <- function(res,
                             TE_feature,
                             plot.title = NULL,
                             output_folder = ".",
                             height = 7,
                             width = 7,
                             device = "png",
                             prefix = NULL) {
  
  # Check for webr package
  if (!requireNamespace("webr", quietly = TRUE)) {
    stop(
      "Package 'webr' is required for donut plots. ",
      "Install it with: ",
      "devtools::install_github('cardiomoon/moonBook')",
      "devtools::install_github('cardiomoon/webr')",
      call. = FALSE
    )
  }

  # Input validation
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res)) {
    stop("'res' must be a data frame.", call. = FALSE)
  }
  
  if (missing(TE_feature)) {
    stop("Argument 'TE_feature' is missing with no default.", call. = FALSE)
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
  
  # Validate device
  if (!device %in% c("png", "pdf", "jpeg", "tiff")) {
    stop(
      "'device' must be one of: png, pdf, jpeg, tiff",
      call. = FALSE
    )
  }
  
  # Get TE_feature type (class, family, or name)  
  TE_type <- .determine_broad_type(res, TE_feature)
  
  if (is.na(TE_type)) {
    stop(
      "TE feature '", TE_feature, "' not found in data.",
      call. = FALSE
    )
  }
  
  # Filter by TE_feature
  res_filtered <- res %>% 
    dplyr::filter(.data[[TE_type]] == TE_feature)
  
  if (nrow(res_filtered) == 0L) {
    stop(
      "No data found for TE feature '", TE_feature, "'.",
      call. = FALSE
    )
  }
  
  # Rename column for webr
  names(res_filtered)[names(res_filtered) == "TE_expression"] <- TE_feature

  # Drop unused factor levels
  res_filtered$expression_type <- droplevels(res_filtered$expression_type)
  
  
  # Create pie plots
  if(!is.null(prefix) & !startsWith(prefix, "_")) {
    prefix <- paste0("_", prefix)
  }
  
  filename <- file.path(output_folder, 
                        paste0("donut_pie_", TE_feature, prefix, ".png"))
  
  # Open device
  if (device == "png") {
    grDevices::png(filename, width = width, height = height,
                   res = 300, units = "in")
  } else if (device == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (device == "jpeg") {
    grDevices::jpeg(filename, width = width, height = height,
                    res = 300, units = "in", quality = 95)
  } else if (device == "tiff") {
    grDevices::tiff(filename, width = width, height = height,
                    res = 300, units = "in")
  }
  
  # Create donut plot
  tryCatch(
    {
      webr::PieDonut(
        res_filtered,
        ggplot2::aes(pies = !!rlang::sym(TE_feature),
                     donuts = !!quote(expression_type)),
        showRatioThreshold = 0.001,
        labelposition = 0,
        pieLabelSize = 8,
        donutLabelSize = 5,
        showPieName = TRUE,
        ratioByGroup = FALSE,
        titlesize = 10
      )
    },
    error = function(e) {
      grDevices::dev.off()
      stop("Failed to create donut plot: ", e$message, call. = FALSE)
    }
  )
  
  grDevices::dev.off()
  
  invisible(NULL)
}
