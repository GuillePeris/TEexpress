#' Create Stacked Bar Plot for Region or TE Class Distribution
#'
#' Internal function to create stacked bar plots showing percentage distribution
#' of genomic regions or TE classes across samples.
#'
#' @param aCR Data frame with columns: region (or class), percent, clon
#' @param plot.title Character string. Plot title.
#' @param analysis Character vector. Labels for x-axis (sample names)
#' @param plot.type Character. Either "region" or "TEclass" to determine
#'   categories and color scheme.
#'
#' @importFrom rlang .data
#' @return A ggplot object.
#'
#' @keywords internal
#' @noRd
stackBarPlot <- function(aCR,
                         plot.title,
                         analysis,
                         plot.type = c("region", "TEclass")) {
  
  plot.type <- match.arg(plot.type)
  
  # Input validation
  if (!is.data.frame(aCR)) {
    stop("'aCR' must be a data frame.", call. = FALSE)
  }
  
  required_cols <- c("region", "percent", "clon")
  if (!all(required_cols %in% colnames(aCR))) {
    stop(
      "'aCR' must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Define levels and colors based on plot type
  if (plot.type == "region") {
    my.levels <- c(
      "Promoter", "5' UTR", "Exon", "Intron", "3' UTR",
      "Downstream", "Intergenic"
    )
    palette <- c(
      "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
      "#EE9B00", "#BB3E03", "#870000"
    )
  } else {  # TEclass
    my.levels <- c("LINE", "LTR", "SINE", "DNA", "Other")
    palette <- c(
      "#e64b35", "#4dbbd5", "#00a087", "#3c5488", "#f39b7f"
    )
    
    # Collapse non-standard classes to "Other"
    aCR$region[!(aCR$region %in% my.levels[1:4])] <- "Other"
  }
  
  names(palette) <- my.levels
  
  # Set factor levels
  aCR$region <- factor(aCR$region, levels = my.levels)
  
  # Create plot using internal helper
  p <- .create_stack_bar(
    data = aCR,
    x = "clon",
    y = "percent",
    fill = "region",
    palette = palette,
    x_labels = analysis,
    title = plot.title
  )
  
  p
}

#' Internal Helper for Creating Stack Bar Plots
#'
#' @keywords internal
#' @noRd
.create_stack_bar <- function(data, x, y, fill, palette, x_labels, title) {
  
  # Create base plot
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(fill = .data[[fill]], y = .data[[y]], x = .data[[x]])
  ) +
    ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.95) +
    ggplot2::scale_x_discrete(
      breaks = levels(data[[x]]),
      labels = x_labels,
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(values = palette)
  
  # Apply theme
  p <- p +
    ggplot2::theme(
      # Grid
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      aspect.ratio = 1,
      
      # Background
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Axis
      axis.text = ggplot2::element_text(colour = "grey30", size = 16),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        size = 20,
        colour = "grey30",
        hjust = 0.5,
        margin = ggplot2::margin(r = 10)
      ),
      axis.ticks = ggplot2::element_line(color = "grey30", linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(0.25, "cm"),
      
      # Legend
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 24, color = "red")
    ) +
    ggplot2::labs(
      title = title,
      y = "% TE loci"
    )
  
  # Format y-axis as percentages
  p <- p +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = c(0, 0)
    )
  
  # Add percentage labels
  p <- p +
    ggfittext::geom_fit_text(
      ggplot2::aes(label = sprintf("%0.1f%%", .data[[y]])),
      stat = "identity",
      min.size = 8,
      size = 16,
      position = ggplot2::position_fill(vjust = 0.5),
      contrast = TRUE,
      show.legend = FALSE
    )
  
  p
}
