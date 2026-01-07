#' Import and Standardize GTF/GFF Annotation File
#'
#' Imports a GTF or GFF annotation file and standardizes it by keeping only
#' standard chromosomes and converting to UCSC chromosome naming style
#' (e.g., "chr1" instead of "1").
#'
#' @param gtf.file Character string. Path to a GTF or GFF format annotation file.
#'   Can be gzip-compressed (.gz).
#' @param format Character string. Format of the annotation file. Must be either
#'   "gtf" or "gff" (also accepts "gff3"). Default is "gtf".
#' @param keep.standard.chroms Logical. If TRUE (default), removes non-standard
#'   chromosomes (keeps only autosomes, sex chromosomes, and mitochondria).
#'   Scaffolds, alternative haplotypes, and patches are removed.
#' @param seqlevels.style Character string. Chromosome naming style to use.
#'   Options include "UCSC" (chr1, chr2, ...), "NCBI" (1, 2, ...), "Ensembl"
#'   (1, 2, ...). Default is "UCSC". Use NULL to skip style conversion.
#'
#' @details
#' This function performs three operations:
#' \enumerate{
#'   \item Imports the GTF/GFF file using \code{rtracklayer::import}
#'   \item Optionally filters to standard chromosomes (autosomes, sex chromosomes,
#'     mitochondria) using \code{GenomeInfoDb::keepStandardChromosomes}
#'   \item Optionally converts chromosome naming to specified style using
#'     \code{GenomeInfoDb::seqlevelsStyle}
#' }
#'
#' @return A \code{GRanges} object containing the genomic annotations with
#'   standardized chromosome names and filtered to standard chromosomes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Import GTF file with default settings
#' gtf <- importGTF("annotations.gtf")
#'
#' # Import GFF3 file
#' gff <- importGTF("annotations.gff3", format = "gff3")
#'
#' # Import without filtering chromosomes
#' gtf_all <- importGTF("annotations.gtf", keep.standard.chroms = FALSE)
#'
#' # Import and keep NCBI chromosome style (1, 2, 3...)
#' gtf_ncbi <- importGTF("annotations.gtf", seqlevels.style = "NCBI")
#'
#' # Import for a specific species
#' gtf_mouse <- importGTF(
#'   "mouse_annotations.gtf",
#'   species = "Mus_musculus"
#' )
#'
#' # Import compressed file
#' gtf <- importGTF("annotations.gtf.gz")
#' }
#'
importGTF <- function(gtf.file,
                      format = "gtf",
                      keep.standard.chroms = TRUE,
                      seqlevels.style = "UCSC") {
  
  # Input validation 
  if (missing(gtf.file)) {
    stop("Argument 'gtf.file' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(gtf.file)) {
    stop("File '", gtf.file, "' not found.", call. = FALSE)
  }
  
  gtf <- tryCatch(
    rtracklayer::import(gtf.file, format = format),
    error = function(e) {
      stop(
        "Failed to import file '", gtf.file, "': ",
        e$message,
        call. = FALSE
      )
    }
  )
  
  # Keep standard chromosomes if requested
  if (keep.standard.chroms) {
    gtf <- GenomeInfoDb::keepStandardChromosomes(gtf, pruning.mode="coarse")
  }
  
  # Convert seqlevels style if requested
  if (!is.null(seqlevels.style)) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- seqlevels.style
  }
  
  gtf
}  
