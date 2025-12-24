#' @import ggplot2
# @import dplyr
#' @import tibble
#' @importFrom grDevices cairo_ps
#' @importFrom utils globalVariables
#' @title Functions to plot DESeq2 results
#'
#' @param res Data frame object with DESeq2 results
#' @param maxpadj P-adjusted value for significant features 
#' @param minlfc Value for dysregulated features
#' @param plot.title Text for graph titles
#' @param device Format for graphs ("pdf", "svg", "eps", "png", "tiff", "jpeg"). 
#'               A vector for several formats can be uses: c("svg", "png") 
#' @param output_folder Folder where graphs will be saved
#' @param width Graph width in inches
#' @param height Graph height in inches
#' @examples
#' \dontrun{
#'   graphTools(res, maxpadj, minlfc, device, 
#'              output, plot.title = plot.title)
#' }
#' 
graphTools <- function(res, maxpadj, minlfc, device = "png", 
                       output_folder = ".", width = 7, height = 7, 
                       plot.title = NULL) {
  # Volcano plot
  filename.volcano <- paste0(output_folder, "/volcanoPlot")
  p.volcano <- volcanoPlot(res, maxpadj, minlfc, plot.title)
  printDevice(p.volcano, filename.volcano, device, width = width, height = width)
  
  # MA-plot
  filename.maplot <- paste0(output_folder, "/maPlot")
  p.maplot <- MAPlot(res, maxpadj, minlfc, plot.title)
  printDevice(p.maplot, filename.maplot, device, width = width, height = width)
}


printDevice <- function(plot, basename, device, width, height) {
  
  # Check devices
  formats <- tolower(gsub("\\.", "", device))
  allowed_devices <- c("pdf", "svg", "eps", "png", "tiff", "jpeg")
  
  stopifnot("Invalid Device(s)" = all(formats %in% allowed_devices))
  
  for (fmt in formats) {
    outfile <- paste0(basename, ".", fmt)
    
    if (fmt == "png") {
      ggsave(outfile, plot, width = width, height = height, 
             dpi = 300, device = "png")
    } else if (fmt == "eps") {
      ggsave(outfile, plot, width = width, height = height, 
             device = cairo_ps, fallback_resolution = 300)
    } else {
      ggsave(outfile, plot, width = width, height = height, device = fmt)
    }
  }
  
  # This is a procedure not returning anything
  invisible()
}

volcanoPlot <- function(res, maxpadj, minlfc, plot.title) {
  # Setting ghost variables to NULL to pass check() 
  log2FoldChange <- log10padj <- diffExpressed <- NULL
  
  # Parameter validation
  if (!is.data.frame(res)) {
    res <- as.data.frame(res)
  }
  if (!all(c("log2FoldChange", "padj") %in% colnames(res))) {
    stop("El data frame debe contener las columnas 'log2FoldChange' y 'padj'")
  }

  # Compute axis limits
  res$log10padj <- -log10(res$padj)
  maxY <- max(res$log10padj[is.finite(res$log10padj)], na.rm = TRUE)
  res$log10padj[!is.finite(res$log10padj)] <- maxY
  maxX <- max(2, max(abs(res$log2FoldChange), na.rm = TRUE) * 1.05)
  maxY <- max(5, maxY*0.90) 
  
  # Assign dysregulated tags
  res$diffExpressed <- dplyr::case_when(
    res$log2FoldChange > minlfc & res$padj < maxpadj ~ "UP",
    res$log2FoldChange < -minlfc & res$padj < maxpadj ~ "DOWN",
    TRUE ~ "NS"
  )
  
  # Counting differentially expressed features
  nDiffExpressed <- c(
    UP = sum(res$diffExpressed == "UP", na.rm = TRUE),
    DOWN = sum(res$diffExpressed == "DOWN", na.rm = TRUE)
  )

  
  # Managing outlier representation
  outPointsList <- res$log10padj > maxY * 0.9
  
  outsidePoints <- NULL
  if (any(outPointsList, na.rm = TRUE)) {
    outsidePoints <- res[outPointsList, ]
    outsidePoints$log10padj <- maxY * 0.9
    res <- res[!outPointsList, ]
  }
  
  # Color management
  mycolors <- c("DOWN" = "#0A9396", "UP" = "#AE2012", "NS" = "grey60")
  
  # Plot base
  p <- ggplot(res, aes(x = log2FoldChange, y = log10padj, 
                       col = diffExpressed)) + 
    geom_point(alpha = 0.7) + 
    theme_classic(base_size = 12) + #, base_family = "Arial") +
    scale_colour_manual(values = mycolors) +
    xlim(-maxX, maxX) + 
    ylim(0, maxY)
  
  # Outliers in y axis
  if (!is.null(outsidePoints)) {
    p <- p + geom_point(data = outsidePoints, 
                        aes(x = log2FoldChange, y = log10padj), 
                        shape = 25, size = 2)
  }
  
  # Reference lines
  p <- p + 
    geom_vline(xintercept = c(-minlfc, minlfc), col = "brown", linetype = 2, linewidth = 0.5) +
    geom_vline(xintercept = 0, linewidth = 0.3)
  
  if (-log10(maxpadj) < maxY) {
    p <- p + geom_hline(yintercept = -log10(maxpadj), col = "brown", linewidth = 0.5)
  }

  # Theme management
  p <- p + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      aspect.ratio = 1, 
      axis.text = element_text(colour = 1, size = 16),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, color = "red")
    ) +
    labs(
      title = plot.title, 
      x = "log2(fold-change)", 
      y = "-log10(adjusted p-value)"
    )
  
  # Text annotations
  p <- p + 
    annotate("text", x = -maxX * 0.98, y = 0.97 * maxY, 
             fontface = 3, size = 4, # family = "Arial",
             label = paste0("Down: ", nDiffExpressed["DOWN"]), 
             colour = "black", hjust = 0) +
    annotate("text", x = maxX * 0.98, y = 0.97 * maxY, 
             fontface = 3, size = 4, # family = "Arial",
             label = paste0("Up: ", nDiffExpressed["UP"]), 
             colour = "black", hjust = 1)
  
  p
}


MAPlot <- function(res, maxpadj, minlfc, plot.title) {
  # Setting ghost variables to NULL to pass check() 
  log2FoldChange <- baseMean <- diffExpressed <- lfc <- NULL
  
  # Parameter validation
  if (!is.data.frame(res)) {
    res <- as.data.frame(res)
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
  # res$diffExpressed <- factor(res$diffExpressed, levels = c("UP", "DOWN", "NS"))
  
  # Color management
  mycolors <- c("DOWN" = "#0A9396", "UP" = "#AE2012", "NS" = "grey60")
  
  # Counting differentially expressed features
  nDiffExpressed <- c(
    UP = sum(res$diffExpressed == "UP", na.rm = TRUE),
    DOWN = sum(res$diffExpressed == "DOWN", na.rm = TRUE)
  )
  
  # Y axis outliers
  outPointsList <- abs(res$log2FoldChange) > maxY
  outsidePoints <- NULL
  if (any(outPointsList, na.rm = TRUE)) {
    outsidePoints <- res[outPointsList, ]
    outsidePoints$lfc <- maxY * sign(outsidePoints$log2FoldChange)
    res <- res[!outPointsList, ]
  }
  
  # Basic plot
  p <- ggplot(res, aes(x = baseMean, y = log2FoldChange, col = diffExpressed)) +
    geom_point(alpha = 0.7) +
    scale_colour_manual(values = mycolors) +
    scale_x_log10() +
    theme_classic(base_size = 12) #, base_family = "Arial")
  
  # Outliers
  if (!is.null(outsidePoints)) {
    p <- p + geom_point(data = outsidePoints, 
                        aes(x = baseMean, y = lfc), 
                        shape = 25, size = 2)
  }
  
  # Reference line
  p <- p + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)
  
  # Y limits
  p <- p + coord_cartesian(ylim = c(-maxY, maxY))
  
  # Theme management
  p <- p + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      aspect.ratio = 1,
      axis.text = element_text(colour = 1, size = 16),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, color = "red")
    ) +
    labs(
      title = plot.title,
      y = "log2(fold-change)",
      x = "mean of normalized counts"
    )
  
  # Annotation position
  max_baseMean <- max(res$baseMean, na.rm = TRUE)
  
  # Annotations
  p <- p + 
    annotate("text", x = max_baseMean, y = -maxY * 0.95, 
             fontface = 3, size = 4, # family = "Arial",
             label = paste0("Underexpressed: ", nDiffExpressed["DOWN"]), 
             colour = "black", hjust = 1) +
    annotate("text", x = max_baseMean, y = maxY * 0.95, 
             fontface = 3, size = 4, # family = "Arial",
             label = paste0("Overexpressed: ", nDiffExpressed["UP"]), 
             colour = "black", hjust = 1)
  p
}