#' Save ggplot to Multiple File Formats
#'
#' Internal function to save a ggplot object to one or more file formats.
#' Supports common vector (SVG, EPS) and raster (PNG, TIFF, JPEG) formats.
#'
#' @param plot A ggplot object to save.
#' @param basename Character string. Base name for output files (without extension).
#'   Extensions will be added automatically based on requested formats.
#' @param device Character vector. File format(s) to save. Supported formats:
#'   "svg", "eps", "png", "tiff", "jpeg". Case-insensitive. Dots are
#'   automatically removed (e.g., ".png" becomes "png").
#' @param width Numeric. Plot width in inches. Default is 7.
#' @param height Numeric. Plot height in inches. Default is 7.
#' 
#'
#' @export
#'


printDevice <- function(plot, basename, device, width = 7, height = 7) {
  
  # Input validation
  if (missing(plot)) {
    stop("Argument 'plot' is missing with no default.", call. = FALSE)
  }
  
  if (missing(basename)) {
    stop("Argument 'basename' is missing with no default.", call. = FALSE)
  }
  
  # Input validation
  if (missing(device)) {
    stop("Argument 'device' is missing with no default.", call. = FALSE)
  }
  
  # Check if output directory exists
  output_dir <- dirname(basename)
  if (output_dir != "." && !dir.exists(output_dir)) {
    stop("Output directory '", output_dir, "' does not exist.", call. = FALSE)
  }
  
  # Check devices
  formats <- tolower(trimws(gsub("\\.", "", device)))  
  allowed_devices <- c("svg", "eps", "png", "tiff", "jpeg")
  
  stopifnot("Invalid Device(s)" = all(formats %in% allowed_devices))
  
  for (fmt in formats) {
    outfile <- paste0(basename, ".", fmt)
    
    tryCatch(
      {
        if (fmt == "png") {
          suppressWarnings(ggplot2::ggsave(
            outfile, plot,
            bg = 'white',
            width = width, height = height,
            dpi = 300, device = "png"
          ))
        } else if (fmt == "eps") {
          suppressWarnings(ggplot2::ggsave(
            outfile, plot,
            bg = 'white',
            width = width, height = height,
            device = grDevices::cairo_ps,
            fallback_resolution = 300
          ))
        } else {
          suppressWarnings(ggplot2::ggsave(
            outfile, plot,
            bg = 'white', 
            width = width, height = height,
            device = fmt
          ))
        }
      },
      error = function(e) {
        warning(
          "Failed to save plot as '", outfile, "': ", e$message,
          call. = FALSE
        )
      }
    )
  }
  
  # This is a procedure not returning anything
  invisible()
}


#' Create Volcano Plot for Differential Expression Results
#'
#' Generates a volcano plot showing log2 fold changes vs. adjusted p-values
#' from differential expression analysis. Points are colored by significance
#' and fold change thresholds.
#'
#' @param res Data frame or DESeqResults object containing differential expression
#'   results. Must include columns: \code{log2FoldChange} and \code{padj}.
#' @param maxpadj Numeric. Adjusted p-value threshold for significance.
#'   Default is 0.05. Features with padj < maxpadj are considered significant.
#' @param minlfc Numeric. Minimum absolute log2 fold change threshold for
#'   significance. Default is 1. Features with |log2FC| > minlfc
#'   are considered differentially expressed.
#' @param plot.title Character string. Title for the plot. Default is
#'   "Volcano Plot".
#'
#' @importFrom rlang .data
#' @return A ggplot object that can be further customized or saved.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Assume 'results' is a DESeq2 results object
#' p <- volcanoPlot(results, maxpadj = 0.05, minlfc = 1,
#'                  plot.title = "Treatment vs Control")
#' print(p)
#'
#' # Save to file
#' printDevice(p, "volcano", "png", width = 8, height = 8)
#' }
volcanoPlot <- function(res,
                        maxpadj = 0.05,
                        minlfc = 1,
                        plot.title = "Volcano Plot") {
  
  
  # Input validation
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }

  if (!is.data.frame(res)) {
    if (inherits(res, "DESeqResults")) {
      res <- as.data.frame(res)
    } else {
      stop("'res' must be a data frame or DESeqResults object.", call. = FALSE)
    }
  }
  
  if (!all(c("log2FoldChange", "padj") %in% colnames(res))) {
    stop("'res' must contain columns: 'log2FoldChange' y 'padj'")
  }

  # Compute -log10(padj)
  res$log10padj <- -log10(res$padj)
  
  # Handle infinite values (padj = 0)
  maxY <- max(res$log10padj[is.finite(res$log10padj)], na.rm = TRUE)
  res$log10padj[!is.finite(res$log10padj)] <- maxY
  
  # Compute axis limits with reasonable scaling
  maxX <- max(2, max(abs(res$log2FoldChange), na.rm = TRUE) * 1.05)
  maxY <- max(5, maxY*0.90) 
  
  # Assign differential expression status  
  res$diffExpressed <- dplyr::case_when(
    res$log2FoldChange > minlfc & res$padj < maxpadj ~ "UP",
    res$log2FoldChange < -minlfc & res$padj < maxpadj ~ "DOWN",
    TRUE ~ "NS"
  )
  
  # Count differentially expressed features
  nDiffExpressed <- c(
    UP = sum(res$diffExpressed == "UP", na.rm = TRUE),
    DOWN = sum(res$diffExpressed == "DOWN", na.rm = TRUE)
  )

  
  # Identify and handle y-axis outliers (top 10%)
  outlier_threshold <- maxY * 0.9
  is_outlier <- res$log10padj > outlier_threshold
  
  outsidePoints <- NULL
  if (any(is_outlier, na.rm = TRUE)) {
    outsidePoints <- res[is_outlier, ]
    outsidePoints$log10padj <- maxY * 0.9
    res <- res[!is_outlier, ]
  }
  
  # Define colors
  color_palette <- c(
    "DOWN" = "#0A9396",  # Teal
    "UP" = "#AE2012",    # Red
    "NS" = "grey60"      # Grey
  )
  
  # Base plot
  p <- ggplot2::ggplot(res, ggplot2::aes(x = .data$log2FoldChange, 
                                         y = .data$log10padj, 
                                         col = .data$diffExpressed)) + 
    ggplot2::geom_point(alpha = 0.7) + 
    ggplot2::theme_classic(base_size = 12) + #, base_family = "Arial") +
    ggplot2::scale_colour_manual(values = color_palette) +
    ggplot2::xlim(-maxX, maxX) + 
    ggplot2::ylim(0, maxY)
  
  # Add outlier points as triangles  
  if (!is.null(outsidePoints)) {
    p <- p + 
      ggplot2::geom_point(
        data = outsidePoints, 
        ggplot2::aes(x = .data$log2FoldChange, y = .data$log10padj), 
        shape = 25, 
        size = 2)
  }
  
  # Reference lines
  p <- p + 
    ggplot2::geom_vline(
      xintercept = c(-minlfc, minlfc), 
      col = "brown", 
      linetype = 2, 
      linewidth = 0.5
      ) +
    ggplot2::geom_vline(
      xintercept = 0, 
      linewidth = 0.3
      )
  
  # Horizontal line for p-value threshold (if visible)
  padj_line_y <- -log10(maxpadj)
  if (padj_line_y < maxY) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = padj_line_y,
        col = "brown",
        linetype = "dashed",
        linewidth = 0.5
      )
  }
  
  # Theme management
  p <- p + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      aspect.ratio = 1, 
      axis.text = ggplot2::element_text(colour = 1, size = 16),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "red")
    ) +
    ggplot2::labs(
      title = plot.title, 
      x = "log2(fold-change)", 
      y = "-log10(adjusted p-value)"
    )
  
  # Text annotations
  p <- p + 
    ggplot2::annotate(
      "text", 
      x = -maxX * 0.98, 
      y = 0.97 * maxY, 
      fontface = 3, 
      size = 4, 
      label = paste0("Down: ", nDiffExpressed["DOWN"]), 
      colour = "black", 
      hjust = 0
      ) +
    ggplot2::annotate(
      "text", 
      x = maxX * 0.98, 
      y = 0.97 * maxY, 
      fontface = 3, 
      size = 4, 
      label = paste0("Up: ", nDiffExpressed["UP"]), 
      colour = "black", 
      hjust = 1
      )
  
  p
}


#' Create MA Plot for Differential Expression Results
#'
#' Generates an MA plot showing log2 fold changes vs. mean normalized counts
#' from differential expression analysis.
#'
#' @param res Data frame or DESeqResults object with columns: baseMean,
#'   log2FoldChange, and padj.
#' @param maxpadj Numeric. Adjusted p-value threshold. Default 0.05.
#' @param minlfc Numeric. Minimum absolute log2 fold change threshold. Default 1.
#' @param plot.title Character string. Plot title.
#' @importFrom rlang .data
#'
#' @return A ggplot object.
#'
#' @keywords internal
#' @noRd
MAPlot <- function(res,
                   maxpadj = 0.05,
                   minlfc = 1,
                   plot.title = "MA Plot") {

  # Input validation
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res)) {
    if (inherits(res, "DESeqResults")) {
      res <- as.data.frame(res)
    } else {
      stop("'res' must be a data frame or DESeqResults object.", call. = FALSE)
    }
  }
 
  if (!all(c("log2FoldChange", "padj", "baseMean") %in% colnames(res))) {
    stop("Data frame must include columns 'baseMean', 'log2FoldChange' and 'padj'")
  }
 
  # Compute axis limits
  maxY <- max(2, max(abs(res$log2FoldChange), na.rm = TRUE) * 0.8)

  # Assign dysregulated tags
  res$diffExpressed <- dplyr::case_when(
    res$log2FoldChange > minlfc & res$padj < maxpadj ~ "UP",
    res$log2FoldChange < -minlfc & res$padj < maxpadj ~ "DOWN",
    TRUE ~ "NS"
  )
  
  # Count differentially expressed features
  nDiffExpressed <- c(
    UP = sum(res$diffExpressed == "UP", na.rm = TRUE),
    DOWN = sum(res$diffExpressed == "DOWN", na.rm = TRUE)
  )
  
  # Y axis outliers
  is_outlier <- abs(res$log2FoldChange) > maxY
  outsidePoints <- NULL
  if (any(is_outlier, na.rm = TRUE)) {
    outsidePoints <- res[is_outlier, ]
    outsidePoints$lfc <- maxY * sign(outsidePoints$log2FoldChange)
    res <- res[!is_outlier, ]
  }
  
  # Define colors
  color_palette <- c(
    "DOWN" = "#0A9396",  # Teal
    "UP" = "#AE2012",    # Red
    "NS" = "grey60"      # Grey
  )
  
  # Base plot
  p <- ggplot2::ggplot(res, 
                       ggplot2::aes(x = .data$baseMean, 
                                    y = .data$log2FoldChange, 
                                    col = .data$diffExpressed)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_colour_manual(values = color_palette) +
    ggplot2::scale_x_log10() 
  
  # Add outlier points as triangles  
  if (!is.null(outsidePoints)) {
    p <- p + 
      ggplot2::geom_point(
        data = outsidePoints, 
        ggplot2::aes(x = .data$baseMean, y = .data$lfc), 
        shape = 25, 
        size = 2
        )
  }
  
  # Reference line
  p <- p + 
    ggplot2::geom_hline(
      yintercept = 0, 
      linetype = 2, 
      linewidth = 0.5
      )
  
  # Y limits
  p <- p + 
    ggplot2::coord_cartesian(
      ylim = c(-maxY, maxY)
      )
  
  # Theme management
  p <- p +
    ggplot2::theme_classic(base_size = 12) + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      aspect.ratio = 1,
      axis.text = ggplot2::element_text(colour = 1, size = 16),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "red")
    ) +
    ggplot2::labs(
      title = plot.title,
      y = "log2(fold-change)",
      x = "mean of normalized counts"
    )
  
  # Text annotations
  max_baseMean <- max(res$baseMean, na.rm = TRUE)
  p <- p + 
    ggplot2::annotate(
      "text", 
      x = max_baseMean, 
      y = -maxY * 0.95, 
      fontface = 3, 
      size = 4, 
      label = paste0("Underexpressed: ", nDiffExpressed["DOWN"]), 
      colour = "black", 
      hjust = 1
      ) +
    ggplot2::annotate(
      "text", 
      x = max_baseMean, 
      y = maxY * 0.95, 
      fontface = 3, 
      size = 4, 
      label = paste0("Overexpressed: ", nDiffExpressed["UP"]), 
      colour = "black", 
      hjust = 1
      )
  p
}
