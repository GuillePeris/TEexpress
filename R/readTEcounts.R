#' Read and Merge TE Count Data by Sample
#'
#' Reads transposable element (TE) count files for multiple samples and merges
#' them into a single count matrix. Count files from TElocal (TEtranscripts)
#' should be tab-separated with two columns: feature names and counts. 
#'
#' @param metadata Data frame with sample metadata, as returned by
#'   \code{\link{read_metadata}}. Must contain columns: File, Sample, Group,
#'   and Condition.
#' @param folder Character string. Path to the directory containing the count
#'   files listed in the first column of \code{metadata}. Can be absolute or
#'   relative path.
#'
#' @details
#' Each count file must:
#' \itemize{
#'   \item Be tab-separated
#'   \item Have a header row
#'   \item Contain exactly 2 columns: feature name and count value
#'   \item Use the same feature names across all samples (order can vary)
#' }
#'
#' The function attempts to efficiently merge count data:
#' \itemize{
#'   \item If all files have identical features in the same order, uses fast
#'     column binding
#'   \item If feature order differs, performs full outer join to align features
#'   \item Missing features in some samples are filled with NA (should not occur
#'     in well-formatted data)
#' }
#'
#' @return A data frame with features as row names and samples as columns.
#'   Column names correspond to the Sample column in \code{metadata}.
#'   Column order matches the row order in \code{metadata}.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Read metadata
#' metadata <- read_metadata("metadata.txt")
#'
#' # Read count files from directory
#' countData <- readTEcounts(metadata, folder = "counts")
#'
#' # Verify dimensions
#' dim(countData)
#' head(countData)
#' }
#'

readTEcounts <- function(metadata, folder) {

  # Input validation
  if (missing(metadata)) {
    stop("Argument 'metadata' is missing with no default.", call. = FALSE)
  }
  
  if (missing(folder)) {
    stop("Argument 'folder' is missing with no default.", call. = FALSE)
  }
  
  if (!dir.exists(folder)) {
    stop("Directory '", folder, "' does not exist.", call. = FALSE)
  }
    
  # Construct file paths
  file_paths <- file.path(folder, metadata[, 1])
  
  # Check all files exist before reading
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop(
      "Missing count file(s):\n  ",
      paste(names(missing_files), ": ", missing_files, collapse = "\n  "),
      call. = FALSE
    )
  }
  
  # Read all count files into a list
  count_list <- lapply(seq_along(file_paths), function(i) {
    counts <- data.table::fread(file_paths[i], sep = "\t", header = TRUE) 
    colnames(counts) <- c("feature", metadata$Sample[i]) 
    counts <- counts %>% tibble::column_to_rownames("feature")
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
  
 count_df
}