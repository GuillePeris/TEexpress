#' @import tibble
#' @title Read TE count data and merge in a data frame by sample
#'
#' @param metadata Data frame with sample information, previously read with
#'        function read_metadata.
#' @param folder Folder for count files listed in first column
#'
#' @returns Data frame with counts by features in rows and samples in columns
#' @export 
#' @examples
#' \dontrun{
#' folder <- "extdata"
#' metadata <- read_metadata(datafile)
#' countData <- readTEcounts(metadata, folder)
#' }
#' 

readTEcounts <- function(metadata, folder) {
  
  # Construct file paths
  file_paths <- file.path(folder, metadata$File)
  
  # Check all files exist before reading
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop("Missing count files:\n  ", paste(missing_files, collapse = "\n  "))
  }
  
  # Read all count files into a list
  count_list <- lapply(seq_along(file_paths), function(i) {
    counts <- data.table::fread(file_paths[i], sep = "\t", header = TRUE) %>% 
               tibble::column_to_rownames(names(.data)[1])
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
  
 count_df
}