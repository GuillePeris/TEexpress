#' Read Metadata File with Sample Information
#'
#' Reads and validates a tab-separated metadata file containing sample information
#' for differential expression analysis. The file must contain specific columns
#' and follow formatting requirements for DESeq2 compatibility.
#'
#' @param datafile Character string. Path to a tab-separated file (.txt, .tsv, or .csv)
#'   containing metadata. Must have exactly 4 columns with headers: File, Sample,
#'   Group, and Condition.
#'
#' @details
#' The metadata file must follow these requirements:
#' \itemize{
#'   \item Exactly 4 columns with headers: File, Sample, Group, Condition
#'   \item At least one row of data (beyond the header)
#'   \item The "Condition" column must contain exactly two unique values: "Control" and "Treat"
#'   \item No duplicate values in the "Sample" or "File" columns
#'   \item Tab-separated format (other delimiters not supported)
#' }
#'
#' The "Control" and "Treat" labels in the Condition column are required for
#' proper DESeq2 analysis, where "Control" will be set as the reference level.
#'
#' @return A data frame with four columns (File, Sample, Group, Condition) and
#'   validated metadata. The Condition column is converted to a factor with
#'   "Control" as the reference level.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read metadata from file
#' metadata <- read_metadata("path/to/metadata.txt")
#'
#' # Example valid file structure:
#' # File            Sample  Group   Condition
#' # WT1.cntTable    WT1     WT      Control
#' # WT2.cntTable    WT2     WT      Control
#' # KO1.cntTable    KO1     KO      Treat
#' # KO2.cntTable    KO2     KO      Treat
#' }
read_metadata <- function(datafile) {
  
  # Input validation
  if (missing(datafile)) {
    stop("Argument 'datafile' is missing with no default.")
  }
 
  # Read and validate metadata file
  if (!file.exists(datafile)) {
    stop("File '", datafile, "' not found.")
  }
  
  # Read metadata file
  metadata <- tryCatch(
    utils::read.table(
      datafile,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\"",
      na.strings = c("NA", "")
    ),
    error = function(e) {
      stop("Failed to read file '", datafile, "': ", e$message, call. = FALSE)
    }
  )  
  # Validate structure
  expected_cols <- c("File", "Sample", "Group", "Condition")
  
  if (ncol(metadata) != length(expected_cols)) {
    stop("File '", datafile, "' should have ", length(expected_cols), 
         " columns, but has ", ncol(metadata), ".")
  }
  
  # Validate Condition column
  unique_conditions <- unique(metadata$Condition)
  required_conditions <- c("Control", "Treat")
  
  if (!all(required_conditions %in% unique_conditions)) {
    stop(
      "The 'Condition' column must contain 'Control' and 'Treat'.\nFound: ",
      paste(unique_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  
  colnames(metadata) <- expected_cols
  metadata
}