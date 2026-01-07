#' Perform Differential Expression Analysis Using DESeq2
#'
#' Performs differential expression analysis on count data using DESeq2.
#' Optionally uses only gene features (non-TE) for size factor estimation.
#'
#' @param countData Data frame with count matrix as returned by
#'   \code{\link{readTEcounts}}. Features in rows, samples in columns.
#'   Row names identify features; TE features are expected to contain ":"
#'   in their names.
#' @param metadata Data frame with sample metadata as returned by
#'   \code{\link{read_metadata}}. Must contain columns: File, Sample, Group,
#'   and Condition. Sample names must match column names in countData.
#' @param useCtrlGenes Logical. If TRUE, only gene features (those without ":"
#'   in their names) are used to estimate size factors. If FALSE, all features
#'   are used. Default is FALSE.
#'
#' @details
#' The function creates a DESeqDataSet with a design formula \code{~ condition},
#' where condition is derived from the Condition column in metadata (Control vs Treat).
#'
#' When \code{useCtrlGenes = TRUE}, size factors are estimated using only
#' gene features (rownames without ":"), which can be more appropriate when
#' analyzing datasets containing both genes and transposable elements. This
#' approach assumes that TEs may have more variable expression and could bias
#' normalization.
#'
#' The function assumes TE features contain ":" in their rownames (e.g.,
#' "transcript_id:TE_name:TE_family:TE_class").
#'
#' @return A DESeqDataSet object after running the DESeq2 analysis pipeline.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read data
#' metadata <- read_metadata("metadata.txt")
#' countData <- readTEcounts(metadata, "counts")
#'
#' # Run DESeq2 using all features for normalization
#' dds <- call_deseq2(countData, metadata, useCtrlGenes = FALSE)
#'
#' # Run DESeq2 using only genes for normalization
#' dds_ctrl <- call_deseq2(countData, metadata, useCtrlGenes = TRUE)
#'
#' # Extract results
#' res <- results_deseq2(dds)
#' }
#'
call_deseq2 <- function(countData, metadata, useCtrlGenes) {
  # Setting ghost variables to NULL to pass check() 
  Condition <- NULL
  
  # Argument validation
  if (missing(countData)) {
    stop("Argument 'countData' is missing with no default.", call. = FALSE)
  }
  
  if (missing(metadata)) {
    stop("Argument 'metadata' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(countData)) {
    stop("'countData' must be a data frame.", call. = FALSE)
  }
  
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame.", call. = FALSE)
  }
  
  
  
  # Check correct factors in Condition.
  expected_levels <- c("Control", "Treat")
  if (any(!metadata$Condition %in% expected_levels)) {
    stop(paste0("Metadata file 'Condition' columns has values ",
                "different to 'Control' and 'Treat'"))
  }
  
  # Create coldata from metadata
  coldata <- metadata[, c("Sample", "Group", "Condition"), drop = FALSE]
  
  # Ensure Condition is a factor with correct levels
  coldata$condition <- factor(metadata$Condition, levels = expected_levels)

  # DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                colData = coldata,
                                design = ~ condition)
  
  # Use genes+TEs or only genes to estimate size factors.
  # We assume TE annotation includes ":"
  # Use genes only or all features for size factor estimation
  if (useCtrlGenes) {
    # Identify gene features (those without ":" in rownames)
    ctrlGenes <- !grepl(":", rownames(dds), fixed = TRUE)
    n_ctrl <- sum(ctrlGenes)
    
    if (n_ctrl == 0) {
      stop("No control genes found (features without ':' in names). ",
           "Cannot use useCtrlGenes = TRUE.", call. = FALSE)
    }
    
    message("Using ", n_ctrl, " gene features (", 
            round(100 * n_ctrl / nrow(dds), 1),
            "%) for size factor estimation...")
    
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = ctrlGenes)
  } else {
    message("Using all ", nrow(dds), " features for size factor estimation...")
  }
  
  # DESeq2 analysis
  dds <- DESeq2::DESeq(dds)
  dds
}