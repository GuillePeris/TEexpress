#' Annotate TEs with Genomic Context Relative to Genes
#'
#' Annotates transposable element (TE) loci with their genomic context relative
#' to protein-coding genes. Uses ChIPseeker to determine if TEs overlap
#' promoters, exons, introns, UTRs, or intergenic regions. Assigns the nearest
#' gene to intergenic TEs.
#'
#' @param df.TEs Data frame containing TE loci with genomic coordinates.
#'   Must include columns: \code{seqnames}, \code{start}, \code{end}, and
#'   \code{strand}. Typically the output from differential expression analysis.
#' @param gene.TxDb A TxDb object containing gene/transcript annotations.
#'   Can be created from GTF files using \code{GenomicFeatures::makeTxDbFromGFF}
#'   or loaded from Bioconductor annotation packages.
#' @param gene_names Data frame mapping gene IDs to gene names. Must contain
#'   columns: \code{gene_id} (matching IDs in TxDb) and \code{gene_name}
#'   (human-readable gene symbols). Used to add gene names to final output.
#' @param transcript.gr GRanges object containing transcript annotations. Must
#'   include a \code{gene_id} column in metadata. Used for assigning the nearest
#'   gene to intergenic TEs. Typically obtained by filtering a gene GTF for
#'   entries where \code{type == "transcript"}.
#' @param TSSminus Integer. Number of bases upstream of TSS to define promoter
#'   region. Should be negative. Default is -5000 (5kb upstream).
#' @param TSSplus Integer. Number of bases downstream of TSS to define promoter
#'   region. Should be positive. Default is 5000 (5kb downstream).
#' @param downstream Integer. Maximum distance downstream of gene end to
#'   consider for "Downstream" annotation. Default is 10000 (10kb).
#'
#' @details
#' The function performs genomic annotation in several steps:
#' \enumerate{
#'   \item Converts TE data frame to GRanges object
#'   \item Annotates TEs using ChIPseeker with priority to coding regions
#'   \item Simplifies annotation categories (e.g., "Exon (uc057wby.1/267097, exon 2 of 7)"
#'     becomes "Exon")
#'   \item For intergenic TEs, assigns the nearest transcript from \code{transcript.gr}
#'   \item Adds gene names from the provided mapping table
#' }
#'
#' Annotation priorities (from highest to lowest):
#' 5' UTR → 3' UTR → Exon → Promoter → Intron → Downstream → Intergenic
#'
#' Note: ChIPseeker may internally convert UCSC chromosome naming (chr1) to
#' NCBI style (1). This function converts them back to UCSC format for consistency.
#'
#'@importFrom rlang .data
#'
#' @return Data frame combining original TE information with genomic annotation
#'   columns:
#' \describe{
#'   \item{annotation}{Simplified genomic context (Promoter, 5' UTR, Exon,
#'     Intron, 3' UTR, Downstream, Intergenic)}
#'   \item{geneChr, geneStart, geneEnd, geneLength, geneStrand}{Genomic
#'     coordinates of the associated gene}
#'   \item{geneId}{Gene identifier from the TxDb}
#'   \item{gene_name}{Human-readable gene symbol (if available in gene_names)}
#'   \item{distanceToTSS}{Distance to transcription start site (from ChIPseeker)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Annotate TEs (assuming TE_results from TE_DEA)
#' TE_annotated <- TE_genomic_context(
#'   df.TEs = TE_results$res.TEs,
#'   gene.TxDb = txdb,
#'   gene_names = gene_names,
#'   transcript.gr = transcript.gr,
#'   TSSminus = -3000,
#'   TSSplus = 3000
#' )
#' }
#'
TE_genomic_context <- function(df.TEs,
                               gene.TxDb,
                               gene_names,
                               transcript.gr,
                               TSSminus = -5000,
                               TSSplus = 5000,
                               downstream = 10000) {

  # Argument validation
  if (missing(df.TEs)) {
    stop("Argument 'df.TEs' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(df.TEs)) {
    stop("'df.TEs' must be a data frame.", call. = FALSE)
  }
  
  required_cols <- c("seqnames", "start", "end", "strand")
  missing_cols <- setdiff(required_cols, colnames(df.TEs))
  
  if (length(missing_cols) > 0L) {
    stop(
      "'df.TEs' must contain columns: ",
      paste(required_cols, collapse = ", "),
      "\nMissing: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (missing(gene.TxDb)) {
    stop("Argument 'gene.TxDb' is missing with no default.", call. = FALSE)
  }
  
  if (!inherits(gene.TxDb, "TxDb")) {
    stop(
      "'gene.TxDb' must be a TxDb object. ",
      "See GenomicFeatures::makeTxDbFromGFF() or use a Bioconductor annotation package.",
      call. = FALSE
    )
  }
  
  # Validate gene_names
  if (missing(gene_names)) {
    stop("Argument 'gene_names' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(gene_names)) {
    stop("'gene_names' must be a data frame.", call. = FALSE)
  }
  
  required_name_cols <- c("gene_id", "gene_name")
  missing_name_cols <- setdiff(required_name_cols, colnames(gene_names))
  
  if (length(missing_name_cols) > 0L) {
    stop(
      "'gene_names' must contain columns: ",
      paste(required_name_cols, collapse = ", "),
      "\nMissing: ", paste(missing_name_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (missing(transcript.gr)) {
    stop("Argument 'transcript.gr' is missing with no default.", call. = FALSE)
  }
  
  if (!inherits(transcript.gr, "GRanges")) {
    stop("'transcript.gr' must be a GRanges object.", call. = FALSE)
  }
  
  # ============================================================
  # Set ChIPseeker Options with Cleanup
  # ============================================================
  options(ChIPseeker.downstreamDistance = downstream)
  options(ChIPseeker.ignore_1st_exon = TRUE)
  options(ChIPseeker.ignore_1st_intron = TRUE)
  
  # ============================================================
  # Convert TEs to GRanges
  # ============================================================
  TE.gr <- tryCatch(
    GenomicRanges::makeGRangesFromDataFrame(df.TEs, keep.extra.columns = FALSE),
    error = function(e) {
      stop(
        "Failed to convert df.TEs to GRanges: ", e$message,
        "\nEnsure columns 'seqnames', 'start', 'end', 'strand' have correct formats.",
        call. = FALSE
      )
    }
  )
  
  # ============================================================
  # Annotate TEs with ChIPseeker
  # ============================================================
  anno <- tryCatch(
    ChIPseeker::annotatePeak(
      TE.gr,
      tssRegion = c(TSSminus, TSSplus),
      TxDb = gene.TxDb,
      # Prioritize coding regions
      genomicAnnotationPriority = c(
        "5UTR", "3UTR", "Exon", "Promoter",
        "Intron", "Downstream", "Intergenic"
      ),
      ignoreDownstream = TRUE,
      # overlap = "all" annotates with actual overlapping gene, not just nearest TSS
      # Also, only consider promoters if overlapping TSS
      overlap = "all"
    ),
    error = function(e) {
      stop("ChIPseeker annotation failed: ", e$message, call. = FALSE)
    }
  )
  
  # Fix Chromosome Naming (ChIPseeker May Convert UCSC to NCBI)  
  anno@anno$geneChr <- ifelse(
    anno@anno$geneChr %in% c("X", "Y", "M", "MT"),
    paste0("chr", anno@anno$geneChr),
    ifelse(grepl("^[0-9]+$", anno@anno$geneChr),
           paste0("chr", anno@anno$geneChr),
           anno@anno$geneChr)
  )
  
  # ============================================================
  # Combine Original Data with Annotations
  # ============================================================
  
  # Preserve row names
  original_rownames <- rownames(df.TEs)
  
  # Add annotation columns
  anno_mcols <- GenomicRanges::mcols(anno@anno)
  df.TEs <- cbind(df.TEs, anno_mcols)
  rownames(df.TEs) <- original_rownames
  df.TEs$distanceToTSS <- NULL
  
  # Simplify annotation labels
  df.TEs <- .simplify_annotations(df.TEs)

  
  # ============================================================
  # Assign Nearest Gene to Intergenic TEs
  # ============================================================

  # Get intergenic TEs
  df.TEs.intergenic <- df.TEs %>%
    dplyr::filter(.data$annotation == "Intergenic") %>%
    dplyr::select(dplyr::all_of(c("seqnames", "start", "end", "strand")))
  
  TEs.intergenic.gr <- tryCatch(
    GenomicRanges::makeGRangesFromDataFrame(df.TEs.intergenic),
    error = function(e) {
      stop(
        "Failed to create GRanges for intergenic TEs: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Find nearest transcript for each intergenic TE
  intergenic_index <- GenomicRanges::nearest(
    TEs.intergenic.gr,
    transcript.gr,
    ignore.strand = FALSE
  )  
  
  # Get nearest transcripts (only for non-NA indices)
  valid_indices <- !is.na(intergenic_index)
  query <- transcript.gr[intergenic_index[valid_indices]]
  
  # Assign gene information to intergenic TEs
  my.columns <- c(
    "geneChr", "geneStart", "geneEnd",
    "geneLength", "geneStrand", "geneId", "transcriptId"
  )
  
  intergenic_rows <- which(df.TEs$annotation == "Intergenic")
  valid_intergenic_rows <- intergenic_rows[valid_indices]
  
  if (length(valid_intergenic_rows) > 0) {
    df.TEs[valid_intergenic_rows, my.columns] <- data.frame(
      geneChr = as.character(GenomicRanges::seqnames(query)),      
      geneStart = BiocGenerics::start(query),
      geneEnd = BiocGenerics::end(query),
      geneLength = BiocGenerics::end(query) - BiocGenerics::start(query),
      geneStrand = as.character(BiocGenerics::strand(query)),
      geneId = GenomicRanges::mcols(query)$gene_id,
      transcriptID = GenomicRanges::mcols(query)$transcript_id
    )
  }
  
  # ============================================================
  # Standardize Strand Notation
  # ============================================================
  
  # ChIPseeker may use numeric strand encoding (1=+, 2=-)
  if ("geneStrand" %in% colnames(df.TEs)) {
    df.TEs$geneStrand <- .convert_strand_notation(df.TEs$geneStrand)
  }

  # ============================================================
  # Add Gene Names
  # ============================================================
  
  # Preserve row names during join
  res.rownames <- rownames(df.TEs)
  
  
  # Check if geneId column exists
  if ("geneId" %in% colnames(df.TEs)) {
    df.TEs <- tryCatch(
      {
        merged <- dplyr::left_join(
          df.TEs,
          gene_names,
          by = c("geneId" = "gene_id")
          # by = dplyr::join_by(geneId == gene_id)
        )
        rownames(merged) <- res.rownames
        merged
      },
      error = function(e) {
        warning(
          "Failed to add gene names: ", e$message,
          "\nReturning data without gene names.",
          call. = FALSE
        )
        df.TEs
      }
    )
  } else {
    warning(
      "'geneId' column not found in annotated data. ",
      "Cannot add gene names.",
      call. = FALSE
    )
  }
  
  df.TEs
}

#' Simplify ChIPseeker Annotation Labels
#'
#' Simplifies detailed ChIPseeker annotations to broad categories.
#' For example, "Exon (uc057wby.1/267097, exon 2 of 7)" becomes "Exon".
#'
#' @param df Data frame with 'annotation' column
#' @return Data frame with simplified annotations
#' @keywords internal
#' @noRd
.simplify_annotations <- function(df) {
  if (!"annotation" %in% colnames(df)) {
    warning("No 'annotation' column found to simplify.", call. = FALSE)
    return(df)
  }
  
  # Create copy to avoid modifying by reference
  df$annotation <- as.character(df$annotation)
  
  # Simplify to broad categories
  # Order matters - more specific patterns first
  df$annotation[startsWith(df$annotation, "Promoter")] <- "Promoter"
  df$annotation[startsWith(df$annotation, "5' UTR")] <- "5' UTR"
  df$annotation[startsWith(df$annotation, "3' UTR")] <- "3' UTR"
  df$annotation[startsWith(df$annotation, "Exon")] <- "Exon"
  df$annotation[startsWith(df$annotation, "Intron")] <- "Intron"
  df$annotation[startsWith(df$annotation, "Downstream")] <- "Downstream"
  df$annotation[startsWith(df$annotation, "Distal")] <- "Intergenic"
  
  df
}

#' Convert Strand Notation to Standard Format
#'
#' Converts numeric strand encoding (1, 2) or other formats to standard (+, -)
#'
#' @param strand Character or numeric vector of strand information
#' @return Character vector with + or - strand notation
#' @keywords internal
#' @noRd
.convert_strand_notation <- function(strand) {
  strand <- as.character(strand)
  strand[!is.na(strand) & strand == "1"] <- "+"
  strand[!is.na(strand) & strand == "2"] <- "-"
  strand
}