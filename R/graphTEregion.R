#' Generate Genomic Region and TE Class Distribution Plots
#'
#' Creates stacked bar plots showing the distribution of transposable elements
#' across genomic regions (promoters, exons, introns, etc.) and TE classes
#' (LINE, SINE, LTR, DNA) for different sample groups and differential
#' expression categories.
#'
#' @param TE_results List object returned by \code{\link{TE_DEA}} and annotated
#'   by \code{\link{annotate_TE_regions}}, containing:
#'   \describe{
#'     \item{res.TEs}{Data frame of DESeq2 results with genomic annotations.
#'       Must include columns: \code{log2FoldChange}, \code{padj},
#'       \code{annotation}, \code{TE_class}}
#'     \item{TE.count}{Data frame of normalized counts with genomic annotations.
#'       Columns should be sample names plus annotation columns}
#'     \item{metadata}{Sample metadata with columns: \code{Sample},
#'       \code{Condition}}
#'   }
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "svg", "eps", "png", "tiff", "jpeg". Can specify multiple formats.
#'   Default is "png".
#' @param output_folder Character string. Directory where plots will be saved.
#'   Default is current directory (".").
#' @param width Numeric. Plot width in inches. Default is 10.
#' @param height Numeric. Plot height in inches. Default is 7.
#' @param minCounts Numeric. Minimum total normalized counts threshold for
#'   considering a TE "expressed" in a condition. TEs below this threshold
#'   across all samples in a condition are excluded. Default is 10.
#' @param minlfc Numeric. Minimum absolute log2 fold change for defining
#'   up/down-regulated TEs. Default is 1.
#' @param maxpadj Numeric. Maximum adjusted p-value for defining significant
#'   differential expression. Default is 0.05.
#' @param plot.title Character string. Title for the plots. If NULL, no title
#'   is added. Default is NULL.
#'
#' @details
#' The function generates two types of stacked bar plots:
#' \enumerate{
#'   \item Genomic Region Distribution: Shows percentage of TEs in each genomic
#'     context (Promoter, 5' UTR, Exon, Intron, 3' UTR, Downstream, Intergenic)
#'   \item TE Class Distribution: Shows percentage of different TE classes
#'     (LINE, SINE, LTR, DNA, Other)
#' }
#'
#' Four categories are plotted for each:
#' \itemize{
#'   \item Control: TEs expressed (total counts > minCounts) in control samples
#'   \item Treatment: TEs expressed in treatment samples
#'   \item Down: Significantly down-regulated TEs (log2FC < -minlfc, padj < maxpadj)
#'   \item Up: Significantly up-regulated TEs (log2FC > minlfc, padj < maxpadj)
#' }
#'
#' Sample group names are derived from the common prefix of sample names in
#' each condition (e.g., "WT1, WT2, WT3" becomes "WT").
#'
#' @section Output Files:
#' Two plot files are created in the specified format(s):
#' \itemize{
#'   \item \code{stackBar_region.<format>}: Genomic region distribution
#'   \item \code{stackBar_TEclass.<format>}: TE class distribution
#' }
#'
#' @importFrom rlang .data
#' @return Invisible NULL. Called for side effects (creating plot files).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After running TE_DEA and TE_regionAnnot
#' TE_results <- TE_DEA(
#'   metafile = "metadata.txt",
#'   folder = "counts",
#'   gtf.TE.file = "TEs.gtf"
#' )
#'
#' TE_results <- TE_regionAnnot(
#'   TE_results = TE_results,
#'   gtf.genes.file = "genes.gtf"
#' )
#'
#' # Generate region distribution plots
#' graphTEregion(
#'   TE_results = TE_results,
#'   device = c("tiff", "png"),
#'   output_folder = "results/plots",
#'   width = 10,
#'   height = 7,
#'   minCounts = 10,
#'   minlfc = 1,
#'   maxpadj = 0.05,
#'   plot.title = "TE Distribution"
#' )
#' }
#' 
graphTEregion <- function(TE_results,
                          device = "png",
                          output_folder = ".",
                          width = 10,
                          height = 7,
                          minCounts = 10,
                          minlfc = 1,
                          maxpadj = 0.05,
                          plot.title = NULL) {

  # Input validation
  if (missing(TE_results)) {
    stop("Argument 'TE_results' is missing with no default.", call. = FALSE)
  }
  
  if (!is.list(TE_results)) {
    stop("'TE_results' must be a list.", call. = FALSE)
  }
  
  required_elements <- c("res.TEs", "TE.count", "metadata")
  missing_elements <- setdiff(required_elements, names(TE_results))
  
  if (length(missing_elements) > 0L) {
    stop(
      "'TE_results' must be a list with elements: ",
      paste(required_elements, collapse = ", "),
      "\nMissing: ", paste(missing_elements, collapse = ", "),
      "\nThis should be the output from TE_DEA().",
      call. = FALSE
    )
  }
  
  # Extract components
  res.TEs <- TE_results$res.TEs
  TE.count <- TE_results$TE.count
  metadata <- TE_results$metadata
  
  # Check output folder exists or can be created
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
  
  # ============================================================
  # Extract Sample Groups
  # ============================================================
  
  # Check for Control and Treat conditions
  unique_conditions <- unique(metadata$Condition)
  
  if (!"Control" %in% unique_conditions) {
    stop(
      "Metadata must contain 'Control' in Condition column. ",
      "Found: ", paste(unique_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!"Treat" %in% unique_conditions) {
    stop(
      "Metadata must contain 'Treat' in Condition column. ",
      "Found: ", paste(unique_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Get control and treatment samples
  ctrlSamples <- metadata$Sample[metadata$Condition == "Control"]
  treatSamples <- metadata$Sample[metadata$Condition == "Treat"]
  
  # Derive group names from sample prefixes
  ctrlName <- .get_common_prefix(ctrlSamples)
  treatName <- .get_common_prefix(treatSamples)
  
  if (nchar(ctrlName) == 0L) {
    ctrlName <- "Control"
  }
  
  if (nchar(treatName) == 0L) {
    treatName <- "Treatment"
  }
    
  expressedTE.ctrl <- .get_expressed(TE.count, ctrlSamples, minCounts)
  expressedTE.treat <- .get_expressed(TE.count, treatSamples, minCounts)
 
  # Extract TEs for each category  
  regions.ctrl <- TE.count[expressedTE.ctrl, ]   
  regions.treat <- TE.count[expressedTE.treat, ]
  
  regions.up <- tryCatch(
    res.TEs %>%
      dplyr::filter(.data$log2FoldChange > minlfc, .data$padj < maxpadj),
    error = function(e) {
      stop("Failed to filter upregulated TEs: ", e$message, call. = FALSE)
    }
  )
  
  regions.down <- tryCatch(
    res.TEs %>%
      dplyr::filter(.data$log2FoldChange < -minlfc, .data$padj < maxpadj),
    error = function(e) {
      stop("Failed to filter downregulated TEs: ", e$message, call. = FALSE)
    }
  )
  
  # Define category names and order
  analysis_categories <- c(ctrlName, treatName, "Down", "Up")

  
  # ============================================================
  # Plot 1: Genomic Region Distribution
  # ============================================================
  region_summary <- tryCatch(
    {
      rbind(
        .summarize_regions(regions.ctrl, ctrlName, "annotation"),
        .summarize_regions(regions.treat, treatName, "annotation"),
        .summarize_regions(regions.down, "Down", "annotation"),
        .summarize_regions(regions.up, "Up", "annotation")
      )
    },
    error = function(e) {
      stop(
        "Failed to combine region summaries: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Set factor levels for proper ordering
  region_summary <- as.data.frame(region_summary)
  region_summary$clon <- factor(region_summary$clon, levels = analysis_categories)
  
  # Create plot
  filename.region <- file.path(output_folder, "stackBar_region")
  
  p.region <- tryCatch(
    stackBarPlot(region_summary, plot.title, analysis_categories, "region"),
    error = function(e) {
      stop(
        "Failed to create region distribution plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot 
  tryCatch(
    printDevice(p.region, filename.region, device,
                width = width, height = height),
    error = function(e) {
      warning(
        "Failed to save region distribution plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # ============================================================
  # Plot 2: TE Class Distribution
  # ============================================================
  
  class_summary <- tryCatch(
    {
      rbind(
        .summarize_regions(regions.ctrl, ctrlName, "TE_class"),
        .summarize_regions(regions.treat, treatName, "TE_class"),
        .summarize_regions(regions.down, "Down", "TE_class"),
        .summarize_regions(regions.up, "Up", "TE_class")
      )
    },
    error = function(e) {
      stop(
        "Failed to combine TE class summaries: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Set factor levels for proper ordering
  class_summary <- as.data.frame(class_summary)
  class_summary$clon <- factor(class_summary$clon, levels = analysis_categories)
  
  # Create plot
  filename.class <- file.path(output_folder, "stackBar_TEclass")
  
  p.class <- tryCatch(
    stackBarPlot(class_summary, plot.title, analysis_categories, "TEclass"),
    error = function(e) {
      stop(
        "Failed to create TE class distribution plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot
  tryCatch(
    printDevice(p.class, filename.class, device,
                width = width, height = height),
    error = function(e) {
      warning(
        "Failed to save TE class distribution plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Don't return anything
  invisible()
}

#' Get Common Prefix from Character Vector
#'
#' Finds the longest common prefix shared by all strings in a vector.
#' Used to derive group names from sample names.
#'
#' @param strings Character vector
#' @return Character string containing the common prefix (may be empty)
#' @keywords internal
#' @noRd
.get_common_prefix <- function(strings) {
  if (length(strings) == 0L) {
    return("")
  }
  
  if (length(strings) == 1L) {
    # For single sample, remove trailing numbers
    return(gsub("[0-9]+$", "", strings[1]))
  }
  
  # Find common prefix
  prefix <- strings[1]
  
  for (i in 2:length(strings)) {
    while (!startsWith(strings[i], prefix)) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
      if (nchar(prefix) == 0) {
        return("")
      }
    }
  }
  
  # Remove trailing non-letter characters
  prefix <- sub("[^A-Za-z]+$", "", prefix)
  
  prefix
}

#' Get expressed TEs
#'
#' Finds TEs with more than 'minCounts' counts.
#'
#' @param TE.count  data frame
#' @param samples samples
#' @param mincounts min Counts
#' @return Names of expressed TEs
#' @keywords internal
#' @noRd
.get_expressed <- function(TE.count, samples, minCounts) {
  tryCatch(
    {
      expressed <- TE.count %>%
        dplyr::filter(
          rowSums(dplyr::across(dplyr::all_of(samples)), na.rm = TRUE) > minCounts
        ) %>%
        rownames()
      expressed
    },
    error = function(e) {
      stop(
        "Failed to filter expressed TEs: ", e$message,
        call. = FALSE
      )
    }
  )
}

#' Summarize TE data frame by variable
#'#'
#' @param regions data frame
#' @param category_name regions
#' @param variable column
#' @importFrom rlang .data
#' @return Summary
#' @keywords internal
#' @noRd
.summarize_regions <- function(regions, category_name, variable) {
  if (nrow(regions) == 0L) {
    # Return empty data frame with correct structure
    return(data.frame(
      region = character(0),
      percent = numeric(0),
      clon = character(0)
    ))
  }
  
  tryCatch(
    {
      summary <- regions %>%
        dplyr::group_by(.data[[variable]]) %>% 
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
          percent = .data$n / sum(.data$n) * 100,
          clon = category_name
        ) %>%
        dplyr::rename(region = .data[[variable]]) %>%   # 
        dplyr::select(.data$region, 
                      .data$percent, 
                      .data$clon)
      
      as.data.frame(summary)
    },
    error = function(e) {
      stop(
        "Failed to summarize regions for ", category_name, ": ",
        e$message,
        call. = FALSE
      )
    }
  )
}