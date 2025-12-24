


#' @import rtracklayer
TE_coords <- function(res.TEs, TE_annot.df) {
  # Setting ghost variables to NULL to pass check() 
  transcript_id <- TE_element <- NULL
  
  # Extract TE info 
  tmp <- do.call(rbind, strsplit(rownames(res.TEs), ":", fixed = TRUE))
  colnames(tmp) <- c("TE_element", "TE_name", "TE_family", "TE_class")
  res.TEs <- cbind(res.TEs, tmp)
  
  # Add coordinates from gtf.TE. Gtf must have column 'transcript_id'
  res.rownames <- rownames(res.TEs)
  res.TEs <- dplyr::left_join(res.TEs, dplyr::select(TE_annot.df, seqnames, 
                                       start, end, strand, transcript_id), 
                  by=dplyr::join_by(TE_element == transcript_id))
  rownames(res.TEs) <- res.rownames
  res.TEs
}
