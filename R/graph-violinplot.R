#' Create violin plot for specific TE list
#'
#' @description
#' Generates a violin plot showing log2 fold changes for a specified list of
#' transposable elements (TEs).
#'
#' @param res.TEs Data frame containing TE differential expression results with
#'   columns: log2FoldChange, padj, TE_element, TE_name, TE_family, TE_class
#' @param TE_list Character vector of TE identifiers to plot. All elements should
#'   belong to the same hierarchical level 
#' @param broad_type Character string specifying the TE classification level
#'   (e.g., "TE_class", "TE_family"). If NULL, automatically detected from TE_list
#' @param minlfc Numeric minimum log2 fold change threshold for significance (default: 1)
#' @param maxpadj Numeric maximum adjusted p-value threshold for significance (default: 0.05)
#' @param min.N Integer minimum number of elements required per group (default: 10)
#' @param width Numeric plot width in inches (default: 7)
#' @param height Numeric plot height in inches (default: 7)
#' @param device Character string specifying output device (default: "png")
#' @param output_folder Character string path to output directory (default: ".")
#' @param plot.title Character string for plot title (default: "Violin plot")
#'
#' @return Invisibly returns NULL. Saves plot to file as side effect.
#'
#' @examples
#' \dontrun{
#' violinPlotByTEList(res.TEs, TE_list = c("LINE", "SINE", "LTR", "DNA"),
#'                    minlfc = 1.5, maxpadj = 0.01)
#' 
#' violinPlotByTEList(res.TEs, TE_list = c("Alu", "L1", "ERV1", "SVA"),
#'                    minlfc = 1.5, maxpadj = 0.01)
#'}
#'
#' @export
violinPlotByTEList <- function(res.TEs,
                               TE_list,
                               broad_type = NULL, 
                               minlfc = 1,
                               maxpadj = 0.05,
                               min.N = 10, 
                               width = 7,
                               height = 7,
                               device = "png",
                               output_folder = ".",
                               plot.title = "Violin plot") {
  
  # Input validation
  if (missing(res.TEs)) {
    stop("Argument 'res.TEs' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res.TEs)) {
    stop("'res.TEs' must be a list.", call. = FALSE)
  }
  
  if (missing(TE_list)) {
    stop("Argument 'TE_list' is missing with no default.", call. = FALSE)
  }

  # Find TE type
  broad_type <- .determine_broad_type(res.TEs, TE_list, broad_type)
  
  # Filter and prepare data for TE type
  filtered_data <- .filter_by_te_list(res.TEs, TE_list, broad_type, min.N)

  # Add expression labels
  filtered_data <- .label_expression(filtered_data, minlfc, maxpadj)
  
  # Create and save plot
  filename <- file.path(output_folder, paste("violinPlot", broad_type, sep = "_"))
  subtitle <- paste0("Expression of specific ", broad_type)
  
  .create_and_save_violin_plot(
    data = filtered_data,
    x = broad_type,
    filename = filename,
    plot.title = plot.title,
    subtitle = subtitle,
    width = width,
    height = height,
    device = device
  )
  
  invisible(NULL)
}

#' Create violin plot for specific TE type
#'
#' @description
#' Generates a violin plot showing the most dysregulated transposable elements
#' within a specific TE type, organized by a finer classification level.
#'
#' @param res.TEs Data frame containing TE differential expression results with
#'   columns: log2FoldChange, padj, TE_element, TE_name, TE_family, TE_class
#' @param TE_type Character string specifying the TE type to analyze 
#'     (example: "LTR", "Alu")
#' @param specific_type Character string for finer classification level:
#'   "TE_family" or "TE_name" (default: "TE_name")
#' @param broad_type Character string for broader classification level:
#'   "TE_class" or "TE_family". If NULL, automatically detected from TE_type 
#'   (default: NULL)
#' @param order Character string specifying expression direction: "up", "down", 
#' or "all" (default: "up")
#' @param minlfc Numeric minimum log2 fold change threshold for significance (default: 1)
#' @param maxpadj Numeric maximum adjusted p-value threshold for significance (default: 0.05)
#' @param nTop Integer number of top dysregulated elements to display (default: 6)
#' @param min.N Integer minimum number of elements required per group (default: 10)
#' @param width Numeric plot width in inches (default: 7)
#' @param height Numeric plot height in inches (default: 7)
#' @param device Character string specifying output device (default: "png")
#' @param output_folder Character string path to output directory (default: ".")
#' @param plot.title Character string for plot title (default: "Violin plot")
#'
#' @return Invisibly returns NULL. Saves plot to file as side effect.
#'
#' @examples
#' \dontrun{
#' violinPlotByTEtype(res.TEs, TE_type = "LINE",
#'                    specific_type = "TE_family", nTop = 10) 
#' 
#' #' violinPlotByTEtype(res.TEs, TE_type = "LINE",
#'                    specific_type = "TE_name", nTop = 10)
#'}
#'
#' @export
violinPlotByTEtype <- function(res.TEs,
                       TE_type,
                       specific_type = "TE_name",
                       broad_type = NULL, 
                       order = "up",
                       minlfc = 1,
                       maxpadj = 0.05,
                       nTop = 6,
                       min.N = 10, 
                       width = 7,
                       height = 7,
                       device = "png",
                       output_folder = ".",
                       plot.title = "Violin plot") {
  
  # Input validation
  if (missing(res.TEs)) {
    stop("Argument 'res.TEs' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res.TEs)) {
    stop("'res.TEs' must be a list.", call. = FALSE)
  }
  
  if (missing(TE_type)) {
    stop("Argument 'TE_type' is missing with no default.", call. = FALSE)
  }
    
  # Find and validate TE type
  broad_type <- .determine_broad_type(res.TEs, TE_type, broad_type)
  .validate_te_types(broad_type, specific_type)
  
  # Filter and prepare data for TE type
  filtered_data <- .filter_by_te_type(res.TEs, TE_type, broad_type, specific_type, min.N)
  
  # Add expression labels
  filtered_data <- .label_expression(filtered_data, minlfc, maxpadj)
  
  # Get top dysregulated elements
  levels.specific.TEs <- .topLevels(filtered_data, specific_type, order, nTop)
  subset.TEs <- filtered_data %>%
    dplyr::filter(.data[[specific_type]] %in% levels.specific.TEs)
  subset.TEs[, specific_type] <- factor(subset.TEs[, specific_type], levels = levels.specific.TEs)
  
  # Create and save plot
  filename <- file.path(output_folder, paste("violinPlot", TE_type, specific_type, sep = "_"))
  subtitle <- paste0("Most ", order, "-dysregulated ", specific_type, " in ", TE_type, " ", broad_type)
  
  .create_and_save_violin_plot(
    data = subset.TEs,
    x = specific_type,
    filename = filename,
    plot.title = plot.title,
    subtitle = subtitle,
    width = width,
    height = height,
    device = device
  )

  invisible(NULL)
}


.topLevels <- function(df, column, order, N) {
  if(order == "all") {
    order <- c("up", "down")
  }
  # We need to extract all levels and then paste them to the upregulated levels
  allLevels <- df %>% 
    dplyr::group_by(.data[[column]]) %>% 
    dplyr::distinct(.data[[column]]) %>% 
    dplyr::pull(.data[[column]])
  upLevels <- df %>% 
    dplyr::group_by(.data[[column]]) %>%  
    dplyr::filter(expression %in% order) %>% 
    dplyr::count() %>% 
    dplyr::arrange(dplyr::desc(.data$n)) %>% 
    as.data.frame() %>%  
    dplyr::pull({{column}})   
  otherLevels <- setdiff(allLevels, upLevels)
  utils::head(c(upLevels, otherLevels), n=N)
}

.find_max_column <- function(df, string) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    NA_character_
  }
  
  counts <- vapply(df, function(x) {
    sum(as.character(x) == string, na.rm = TRUE)
  }, integer(1))
  
  if (max(counts) == 0L) {
    return(NA_character_)
  }
  
  names(which.max(counts))
}

#' Determine broad type from TE identifier(s)
#' @noRd
.determine_broad_type <- function(res.TEs, te_identifier, broad_type = NULL) {
  if (is.null(broad_type)) {
    broad_type <- .find_max_column(res.TEs, te_identifier[1])
  }
  
  if (is.na(broad_type)) {
    stop("Element '", te_identifier[1], "' not found", call. = FALSE)
  }
  
  # For lists, check all elements belong to same type
  if (length(te_identifier) > 1) {
    all_types <- vapply(te_identifier, function(te) .find_max_column(res.TEs, te),
                        character(1), USE.NAMES = FALSE)
    
    if(any(is.na(all_types))) {
      stop(
        "Some TEs were not found in results: ",
        paste(te_identifier[is.na(all_types)], collapse = ", "),
        call. = FALSE
      )
    }
    
    if (length(unique(all_types)) > 1L) {
      stop(
        "TEs do not belong to same type. Found types: ",
        paste(unique(all_types), collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  broad_type
}

#' Filter data by TE list
#' @noRd
.filter_by_te_list <- function(res.TEs, TE_list, broad_type, min.N) {
  res.TEs %>%
    dplyr::select(.data$log2FoldChange,
                  .data$padj,
                  .data$TE_element,
                  .data$TE_name,
                  .data$TE_family,
                  .data$TE_class) %>%
    dplyr::filter(.data[[broad_type]] %in% TE_list) %>%
    dplyr::group_by(.data[[broad_type]]) %>%
    dplyr::filter(dplyr::n() >= min.N) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

#' Filter data by TE type
#' @noRd
.filter_by_te_type <- function(res.TEs, TE_type, broad_type, specific_type, min.N) {
  res.TEs %>%
    dplyr::select(.data$log2FoldChange,
                  .data$padj,
                  .data$TE_element,
                  .data$TE_name,
                  .data$TE_family,
                  .data$TE_class) %>%
    dplyr::filter(.data[[broad_type]] == TE_type) %>%
    dplyr::group_by(.data[[specific_type]]) %>%
    dplyr::filter(dplyr::n() >= min.N) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

#' Add expression labels to data
#' @noRd
.label_expression <- function(data, minlfc, maxpadj) {
  data$expression <- "ns"
  data$expression[data$log2FoldChange >= minlfc & data$padj <= maxpadj] <- "up"
  data$expression[data$log2FoldChange <= -minlfc & data$padj <= maxpadj] <- "down"
  data
}

#' Create and save violin plot
#' @noRd
.create_and_save_violin_plot <- function(data, x, filename, plot.title,
                                         subtitle, width, height, device) {
  ylim <- c(min(data$log2FoldChange) - 0.25, max(data$log2FoldChange) + 0.25)
  
  p <- ggpubr::ggviolin(data, x = x, y = "log2FoldChange", alpha = 0.5,
                        fill = "#4d6600", color = "#333300",
                        draw_quantiles = 0.5,
                        short.panel.labs = FALSE, width = 0.6, ylim = ylim,
                        outlier.shape = NA)
  
  p <- p +
    gghalves::geom_half_point(
      side = "l",
      data = data[data$expression == "up", ],
      shape = 21, range_scale = .4,
      alpha = 0.5,
      fill = "#990f02", color = "black", stroke = 0.2,
      show.legend = TRUE, width = 0.25
    ) +
    gghalves::geom_half_point(
      data = data[data$expression == "down", ],
      shape = 21, alpha = 0.5,
      fill = "#0f4392", color = "black", stroke = 0.2,
      show.legend = FALSE, width = 0.25
    )
  
  p <- p +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        face = "bold", color = "#993333",
        size = 12, angle = 45, hjust = 1
      )
    ) +
    ggplot2::labs(y = "log2(FC)") +
    ggplot2::guides(shape = "none") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2) +
    ggplot2::ggtitle(plot.title, subtitle = subtitle) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_y_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6))
  
  # Save plot
  tryCatch(
    printDevice(p, filename, device, width = width, height = height),
    error = function(e) {
      warning(
        "Failed to save violin plot: ", e$message,
        call. = FALSE
      )
    }
  )
}

#' Validate TE type hierarchy
#' @noRd
.validate_te_types <- function(broad_type, specific_type) {
  TE_elements_broad <- c("TE_class", "TE_family")
  if (!(broad_type %in% TE_elements_broad)) {
    stop("'broad_type' must be in ",
         paste(TE_elements_broad, collapse = ", "),
         call. = FALSE)
  }
  
  TE_elements_specific <- c("TE_family", "TE_name")
  if (!(specific_type %in% TE_elements_specific)) {
    stop("'specific_type' must be in ",
         paste(TE_elements_specific, collapse = ", "),
         call. = FALSE)
  }
  
  if (broad_type == specific_type) {
    stop("'specific_type' has to be hierarchically lower than 'broad_type'", call. = FALSE)
  }
}
