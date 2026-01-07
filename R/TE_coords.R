#' Add Genomic Coordinates to TE Count Data
#'
#' Extracts transposable element (TE) information from count data rownames
#' and adds genomic coordinates by joining with TE annotation data.
#'
#' @param countData.TEs Data frame with TE count data. Rownames must follow
#'   the format: "TE_element:TE_name:TE_family:TE_class" (colon-separated).
#'   Each column should contain count data for a sample.
#' @param TE_annot.df Data frame with TE annotations, typically derived from
#'   a GTF file. Must contain columns: seqnames, start, end, strand, and
#'   transcript_id. The transcript_id column is used to match with the first
#'   field in countData.TEs rownames.
#'
#'
#' @return Data frame containing the original count data plus additional columns:
#'   \itemize{
#'     \item TE_element: Transcript ID extracted from rownames
#'     \item TE_name: TE name extracted from rownames
#'     \item TE_family: TE family extracted from rownames
#'     \item TE_class: TE class extracted from rownames
#'     \item seqnames: Chromosome/sequence name from annotation
#'     \item start: Start position from annotation
#'     \item end: End position from annotation
#'     \item strand: Strand information from annotation
#'   }
#'   Original rownames are preserved.
#'   
#'   @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming countData.TEs has rownames like:
#' # "L1PA2_dup501:L1PA2:L1:LINE"
#' # "AluY_dup68349:AluY:Alu:SINE"
#'
#' # And TE_annot.df has columns: transcript_id, seqnames, start, end, strand
#'
#' countData_with_coords <- TE_coords(countData.TEs, TE_annot.df)
#'
#' # Check results
#' head(countData_with_coords)
#' }
#'
TE_coords <- function(countData.TEs, TE_annot.df) {

  # Argument validation
  if (missing(countData.TEs)) {
    stop("Argument 'countData.TEs' is missing with no default.", call. = FALSE)
  }
  
  if (missing(TE_annot.df)) {
    stop("Argument 'TE_annot.df' is missing with no default.", call. = FALSE)
  }  
  
  # Check required columns in TE_annot.df
  required_cols <- c("seqnames", "start", "end", "strand", "transcript_id")
  missing_cols <- setdiff(required_cols, colnames(TE_annot.df))
  if (length(missing_cols) > 0) {
    stop("'TE_annot.df' is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  # Check rownames format (should contain colons)
  res.rownames <- rownames(countData.TEs)
  if (!any(grepl(":", res.rownames, fixed = TRUE))) {
    stop("'countData.TEs' rownames must be in format: ",
         "'transcript_id:TE_name:TE_family:TE_class' (colon-separated).",
         call. = FALSE)
  }
  
  # Split rownames
  tmp <- tryCatch(
    do.call(rbind, strsplit(res.rownames, ":", fixed = TRUE)),
    error = function(e) {
      stop("Failed to split rownames by ':'. ", e$message, call. = FALSE)
    }
  )
  colnames(tmp) <- c("TE_element", "TE_name", "TE_family", "TE_class")
  countData.TEs <- cbind(countData.TEs, tmp)
  
  # Add coordinates from gtf.TE. Gtf must have column 'transcript_id'
  countData.TEs <- dplyr::left_join(countData.TEs, dplyr::select(TE_annot.df, 
                                                                 .data$seqnames, 
                                                                 .data$start, 
                                                                 .data$end, 
                                                                 .data$strand, 
                                                                 .data$transcript_id), 
                  by= c("TE_element" = "transcript_id") ) #dplyr::join_by(TE_element == transcript_id))
  rownames(countData.TEs) <- res.rownames
  countData.TEs
}

#' @importFrom rlang .data
addTEColumns <- function(df, TE.coordinates) {
  df.rownames <- rownames(df)
  tmp <- as.data.frame(do.call(rbind, strsplit(rownames(df), ":", fixed = TRUE)))
  colnames(tmp) <- c("TE_element", "TE_name", "TE_family", "TE_class")
  
  df <- cbind(df, tmp)
  df <- dplyr::left_join(df, dplyr::select(TE.coordinates, 
                                           .data$seqnames, 
                                           .data$start, 
                                           .data$end, 
                                           .data$strand, 
                                           .data$TE_element), 
                         by="TE_element")
  rownames(df) <- df.rownames
  df
}

