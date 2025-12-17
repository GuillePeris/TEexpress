#' @importFrom utils read.table
#' @title Reads metadata file with sample information
#'
#' @param datafile A tab separated file with four columns and a header. Structure 
#' must follow this example:
#' 
#' File         Sample	Group	 Condition
#' WT1.cntTable	WT1	    WT	   Control
#' ....
#' KO1.cntTable	KO1	    KO	   Treat
#' 
#' "Control" and "Treat" labels are required for DESeq2 analysis.
#' 
#' @returns Data frame with metadata
#' @export
#' @examples 
#' \dontrun{
#' datafile <- "extdata/data.csv"
#' metadata <- read_metadata(datafile)
#' }
read_metadata <- function(datafile) {
  # Define expected columns
  expected_cols <- c("File", "Sample", "Group", "Condition")
  
  # Read and validate metadata file
  if (!file.exists(datafile)) {
    stop("File '", datafile, "' not found.")
  }
  
  metadata <- read.table(datafile, sep = "\t", header = TRUE, 
                         check.names = FALSE, stringsAsFactors = FALSE)
  
  # Validate structure
  if (ncol(metadata) != length(expected_cols)) {
    stop("File '", datafile, "' should have ", length(expected_cols), 
         " columns, but has ", ncol(metadata), ".")
  }
  
  colnames(metadata) <- expected_cols
  metadata
}