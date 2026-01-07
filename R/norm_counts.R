#' Get Normalized Counts from DESeq2 Object
#'
#' Extracts normalized count data from a DESeqDataSet object. Counts are
#' normalized by size factors to account for differences in sequencing depth
#' between samples.
#'
#' @param dds A DESeqDataSet object as returned by \code{\link{call_deseq2}}.
#'   Size factors must have been calculated (automatically done by
#'   \code{DESeq2::DESeq()}).
#'
#' @details
#' This is a convenience wrapper around \code{DESeq2::counts(dds, normalized = TRUE)}.
#' 
#' Normalized counts are calculated by dividing raw counts by size factors.
#' These counts are appropriate for visualization and clustering, but should
#' NOT be used as input for differential expression analysis (use raw counts
#' for that).
#'
#' @return A numeric matrix with normalized counts. Features in rows, samples
#'   in columns. Row and column names are preserved from the DESeqDataSet.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run DESeq2 analysis
#' dds <- call_deseq2(countData, metadata)
#'
#' # Extract normalized counts
#' norm_counts_matrix <- norm_counts(dds)
#'
#' # Convert to data frame if needed
#' norm_counts_df <- as.data.frame(norm_counts_matrix)
#' }
norm_counts <- function(dds) {
  if (missing(dds)) {
    stop("Argument 'dds' is missing with no default.", call. = FALSE)
  }
  
  if (!methods::is(dds, "DESeqDataSet")) {
    stop("'dds' must be a DESeqDataSet object.", call. = FALSE)
  }
  
  DESeq2::counts(dds, normalized=TRUE)
}