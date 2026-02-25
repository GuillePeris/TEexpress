#' Get Results After DESeq2 Analysis
#'
#' Extracts differential expression results from a DESeqDataSet object,
#' with optional log2 fold change shrinkage for improved effect size estimates.
#'
#' @param dds A DESeqDataSet object as returned by \code{\link{call_deseq2}} or
#'   \code{DESeq2::DESeq()}. Must have completed the DESeq2 analysis pipeline.
#' @param shrinklog2FC Logical. If TRUE, applies log2 fold change shrinkage using
#'   the apeglm method to produce more accurate effect size estimates. If FALSE,
#'   returns standard DESeq2 results. Default is FALSE. Note: Shrinkage requires
#'   the \code{apeglm} package to be installed.
#'
#' @details
#' This function extracts differential expression results comparing the "Treat"
#' condition against the "Control" condition (reference level).
#'
#' When \code{shrinklog2FC = TRUE}, log2 fold changes are shrunk using the
#' apeglm method (\code{DESeq2::lfcShrink})
#'
#' Rows with NA values (typically low-count genes filtered by DESeq2's
#' independent filtering) are removed, and results are sorted by adjusted p-value.
#'
#' @return A data frame with DESeq2 results with rows sorted by adjusted p-value 
#' (most significant first). Row names correspond to feature names from the count data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run DESeq2 analysis
#' dds <- call_deseq2(countData, metadata, useCtrlGenes = FALSE)
#'
#' # Get results without shrinkage
#' res <- results_deseq2(dds, shrinklog2FC = FALSE)
#'
#' # Get results with log2FC shrinkage (requires apeglm package)
#' res_shrink <- results_deseq2(dds, shrinklog2FC = TRUE)
#'
#' # View top results
#' head(res)
#' }
#'
results_deseq2 <- function(dds, shrinklog2FC=FALSE) {
  # Setting ghost variables to NULL to pass check() 
  colData <- NULL
  
  # Essential validation
  if (missing(dds)) {
    stop("Argument 'dds' is missing with no default.", call. = FALSE)
  }
  
  # Checing levels
  expected_levels <- levels(SummarizedExperiment::colData(dds)$condition)
  if (length(expected_levels) != 2) {
    stop("'condition' must have exactly 2 levels, found ",
         length(expected_levels), ".", call. = FALSE)
  }
  
  # Shrink log2FC if requires
  if(shrinklog2FC) {
    # Check if apeglm is installed
    if (!requireNamespace("apeglm", quietly = TRUE)) {
      stop("Package 'apeglm' is required for log2FC shrinkage but is not installed.\n",
           "Install it with: BiocManager::install('apeglm')\n",
           "Or set shrinklog2FC = FALSE.", call. = FALSE)
    }
    res <- DESeq2::lfcShrink(dds, coef=paste0("condition_", expected_levels[2], 
                                              "_vs_", expected_levels[1]),  
                             type="apeglm")
  } else {
    res <- DESeq2::results(dds, 
                           contrast = c("condition", expected_levels[2], expected_levels[1]))
  } 
  
  res <- as.data.frame(res)
  
  # Remove NA rows and sort by padj
  res_complete <- res[stats::complete.cases(res), ]
  res_sorted <- res_complete[order(res_complete$padj), ]
  res_sorted
}