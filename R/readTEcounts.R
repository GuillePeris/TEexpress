#' Read TE count data and merge in a data frame by sample
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
#' @param folder Folder for count files listed in first column
#'
#' @returns Data frame with counts by features in rows and samples in columns
#' @export 
#'
#' @examples 
#' datafile <- "my.path/data.csv"
#' folder <- "data/"
#' df <- readTEexpress(datafile, folder)

readTEcounts <- function(datafile, folder) {
  
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
  
  # Construct file paths
  file_paths <- file.path(folder, metadata$File)
  
  # Check all files exist before reading
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop("Missing count files:\n  ", paste(missing_files, collapse = "\n  "))
  }
  
  # Read all count files into a list
  count_list <- lapply(seq_along(file_paths), function(i) {
    counts <- read.table(file_paths[i], sep = "\t", header = TRUE, 
                         check.names = FALSE, row.names = 1)
    colnames(counts) <- metadata$Sample[i]
    counts
  })
  
  # Check if all files have identical row names (common case)
  all_rownames <- lapply(count_list, rownames)
  identical_order <- all(sapply(all_rownames[-1], identical, all_rownames[[1]]))
  
  # Merge efficiently based on row name consistency
  if (identical_order) {
    # Fast path: simple column binding
    count_df <- do.call(cbind, count_list)
  } else {
    # Slow path: merge with row name matching
    count_df <- count_list[[1]]
    for (i in seq(2, length(count_list))) {
      count_df <- merge(count_df, count_list[[i]], by = "row.names", all = TRUE)
      rownames(count_df) <- count_df$Row.names
      count_df$Row.names <- NULL
    }
  }
  
  return(count_df)
}