#' Annotate TEs with Genomic Regions and Generate Region Distribution Plots
#'
#' Annotates transposable elements (TEs) with their genomic context relative to
#' protein-coding genes and generates stacked bar plots showing the distribution
#' of TEs across genomic regions (promoters, exons, introns, etc.). This is a
#' wrapper function that combines genomic annotation with visualization.
#'
#' @param TE_results List object returned by \code{\link{TE_DEA}}, containing:
#'   \describe{
#'     \item{res.TEs}{Data frame of DESeq2 results for TEs}
#'     \item{TE.count}{Data frame of normalized counts for TEs}
#'     \item{gene.count}{Data frame of normalized counts for genes}
#'     \item{metadata}{Sample metadata data frame}
#'   }
#' @param gtf.genes.file Character string. Full path to GTF file containing
#'   gene annotations. Should include protein-coding genes with transcript
#'   information. Must have columns: \code{type}, \code{transcript_biotype},
#'   \code{gene_biotype}, \code{gene_id}, \code{gene_name}.
#' @param gtf.format Character string. Format of gene GTF file: "gtf" (default),
#'        "gff3· or "gff".
#' @param output_folder Character string. Path for output directory. A
#'   subdirectory "TEs_annotated" will be created inside this folder. Default
#'   is current directory (".").
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "svg", "eps", "png", "tiff", "jpeg". Can specify multiple formats.
#'   Default is "png".
#' @param plot.title Character string. Title for the region distribution plots.
#'   If NULL (default), a generic title will be used.
#' @param maxpadj Numeric. Adjusted p-value threshold for significance.
#'   Features with padj < maxpadj are considered significant. Default is 0.05.
#' @param minlfc Numeric. Minimum absolute log2 fold change threshold for
#'   differential expression. Default is 1.
#' @param width Numeric. Plot width in inches. Default is 10.
#' @param height Numeric. Plot height in inches. Default is 7.
#' @param minCounts Numeric. Minimum total normalized counts threshold for 
#'   considering a TE loci as expressed. Used for region distribution plots. 
#'   Default is 10.
#' @param is_ext_3UTR Boolean. TRUE if extended 3'UTR analysis is to be
#'        performed. Defaults to FALSE. To be implemented
#' @param ext_3UTR_file Character string. Filename for 3'UTR analysis result,
#'                   in gff format. To be implemented.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Imports and filters gene GTF to protein-coding genes only
#'   \item Creates a TxDb object for ChIPseeker annotation
#'   \item Extracts transcript features for nearest gene assignment
#'   \item Annotates TEs with genomic context using \code{\link{TE_genomic_context}}
#'   \item Saves annotated results to TSV file
#'   \item Generates region distribution plots using \code{\link{graphTEregion}}
#' }
#'
#' Only protein-coding genes are considered for annotation to focus on
#' biologically relevant gene-TE associations. Non-coding genes and pseudogenes
#' are excluded.
#'
#' @section Output Files:
#' The function creates a subdirectory structure:
#' \preformatted{
#' output_folder/
#'   └── TEs_annotated/
#'       ├── DESeq2_TE_results_annotated.tsv
#'       └── [region distribution plots in specified format(s)]
#' }
#' 
#' @importFrom rlang .data
#'
#' @return A list with three updated elements for downstream analysis:
#' \describe{
#'   \item{res.TEs}{Data frame of DESeq2 results with added genomic annotation
#'     columns (annotation, geneChr, geneStart, geneEnd, geneId, gene_name, etc.)}
#'   \item{TE.count}{Data frame of normalized TE counts with added genomic annotation
#'     columns}
#'   \item{gene.count}{Data frame of normalized gene counts}
#'   \item{metadata}{Original sample metadata (unchanged)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run differential expression analysis first
#' TE_results <- TE_DEA(
#'   metafile = "metadata.txt",
#'   folder = "counts",
#'   gtf.TE.file = "TEs.gtf"
#' )
#'
#' # Annotate TEs and generate region plots
#' TE_results_annotated <- annotate_TE_regions(
#'   TE_results = TE_results,
#'   gtf.genes.file = "genes.gtf",
#'   output_folder = "results",
#'   device = c("jpeg", "png"),
#'   plot.title = "TE Distribution by Genomic Region",
#'   minCounts = 10
#' )
#'
#' }
#'
annotate_TE_regions <- function(TE_results,
                           gtf.genes.file,
                           gtf.format = "gtf",
                           output_folder = ".",
                           device = "png",
                           plot.title = NULL,
                           maxpadj = 0.05,
                           minlfc = 1,
                           width = 10,
                           height = 7,
                           minCounts = 10,
                           is_ext_3UTR = FALSE,
                           ext_3UTR_file = NULL) {
  
  message("============================================") 
  message("     TE loci genomic region annotation      ")
  message("============================================") 
  
  # Input validation
  if (missing(TE_results)) {
    stop("Argument 'TE_results' is missing with no default.", call. = FALSE)
  }
  
  if (!is.list(TE_results)) {
    stop("'TE_results' must be a list.", call. = FALSE)
  }
  
  required_elements <- c("res.TEs", "TE.count", "gene.count", "metadata")
  missing_elements <- setdiff(required_elements, names(TE_results))
  
  if (length(missing_elements) > 0L) {
    stop(
      "'TE_results' must be a list with elements: ",
      paste(required_elements, collapse = ", "),
      "\nMissing: ", paste(missing_elements, collapse = ", "),
      "\nThis should be the output from TE_DEA().",
      call. = FALSE
    )
  }
  
  # Validate gtf.genes.file
  if (missing(gtf.genes.file)) {
    stop("Argument 'gtf.genes.file' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(gtf.genes.file)) {
    stop("Gene GTF file '", gtf.genes.file, "' not found.", call. = FALSE)
  }
  
  # Extract and validate components
  res.TEs <- TE_results$res.TEs
  TE.count <- TE_results$TE.count
  gene.count <- TE_results$gene.count
  metadata <- TE_results$metadata
  
  # ============================================================
  # Step 1: Import and Filter Gene GTF
  # ============================================================
  
  message("==> Importing gene GTF file")
  start.time <- Sys.time()
  
  gtf.genes <- tryCatch(
    importGTF(gtf.genes.file, format = gtf.format),
    error = function(e) {
      stop("Failed to import gene GTF file: ", e$message, call. = FALSE)
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished importing GTF (", round(duration[[1]], 2), " seconds)")
  
  # ============================================================
  # Step 2: Filter to Protein-Coding Genes
  # ============================================================
  
  message("==> Filtering to protein-coding genes")
  start.time <- Sys.time()
  
  # Check required columns exist
  gtf_mcols <- GenomicRanges::mcols(gtf.genes)
  required_gtf_cols <- c("type", "transcript_biotype")
  missing_gtf_cols <- setdiff(required_gtf_cols, colnames(gtf_mcols))
  
  if (length(missing_gtf_cols) > 0L) {
    stop(
      "Gene GTF is missing required column(s): ",
      paste(missing_gtf_cols, collapse = ", "),
      "\nRequired columns: ", paste(required_gtf_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Filter to protein-coding genes
  is_protein_coding_gene <- !is.na(gtf.genes$gene_biotype) &
    gtf.genes$gene_biotype == "protein_coding"
  gtf.genes <- gtf.genes[is_protein_coding_gene]
  
  # Filter to protein-coding transcripts
  is_protein_coding <- !is.na(gtf.genes$transcript_biotype) &
    gtf.genes$transcript_biotype == "protein_coding"
  
  if (!any(is_protein_coding)) {
    stop(
      "No protein-coding transcripts found in GTF. ",
      "Check that the GTF has 'transcript_biotype' column with ",
      "'protein_coding' values.",
      call. = FALSE
    )
  }
  
  transcript.gr <- gtf.genes[is_protein_coding]   
  
  # Extract transcript features
  is_transcript <- transcript.gr$type == "transcript"
  
  if (!any(is_transcript)) {
    stop(
      "No transcript features found in filtered GTF. ",
      "The GTF must contain entries with type='transcript'.",
      call. = FALSE
    )
  }
  
  transcript.gr <- transcript.gr[is_transcript]  
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished filtering (", round(duration[[1]], 2), " seconds)")

  # ============================================================
  # Step 3: Create TxDb Object
  # ============================================================
  
  message("==> Creating TxDb object for ChIPseeker")
  start.time <- Sys.time()
  
  gene.TxDb <- tryCatch(
    suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gtf.genes)),   
    error = function(e) {
      stop(
        "Failed to create TxDb from GTF: ", e$message,
        "\nEnsure GTF has proper structure with transcript and exon features.",
        call. = FALSE
      )
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished creating TxDb (", round(duration[[1]], 2), " seconds)")  
  
  # ============================================================
  # Step 4: Prepare Gene Names Table
  # ============================================================

  message("==> Preparing gene name mapping")
  start.time <- Sys.time()
  
  # Check for required columns
  if (!"gene_id" %in% colnames(gtf_mcols) |
      !"gene_name" %in% colnames(gtf_mcols))  {
    stop("GTF must contain 'gene_id' and 'gene_name' columns.", call. = FALSE)
  } else {
    # Filter and prepare gene names
    gene_names <- tryCatch(
      {
        df <- as.data.frame(gtf.genes)
        df %>%
          dplyr::filter(
            .data$type == "gene",
            !is.na(.data$gene_biotype),
            .data$gene_biotype == "protein_coding"
          ) %>%
          dplyr::select(.data$gene_id, .data$gene_name) %>%
          unique()
      },
      error = function(e) {
        stop(
          "Failed to extract gene names: ", e$message,
          "\nProceeding without gene names.",
          call. = FALSE
        )
      }
    )
  }
    end.time <- Sys.time()
    duration <- difftime(end.time, start.time, units = "secs")
    message("       -> Finished preparing gene names (", round(duration[[1]], 2), " seconds)")
  
  # ============================================================
  # Step 5: Annotate TEs with Genomic Context
  # ============================================================
  
  message("==> Annotating TEs with genomic context")
  start.time <- Sys.time()
  
  res.TEs <- tryCatch(
    TE_genomic_context(res.TEs, 
                       gene.TxDb, 
                       gene_names, 
                       transcript.gr,
                       is_ext_3UTR = is_ext_3UTR,
                       ext_3UTR_file = ext_3UTR_file),
    error = function(e) {
      stop(
        "Failed to annotate res.TEs with genomic context: ",
        e$message,
        call. = FALSE
      )
    }
  )
  
  TE.count <- tryCatch(
    TE_genomic_context(TE.count, 
                       gene.TxDb, 
                       gene_names, 
                       transcript.gr,
                       is_ext_3UTR = is_ext_3UTR,
                       ext_3UTR_file = ext_3UTR_file),
    error = function(e) {
      stop(
        "Failed to annotate TE.count with genomic context: ",
        e$message,
        call. = FALSE
      )
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished annotation (", round(duration[[1]], 2), " seconds)")
  
  #============================================================
  # Step 6: Save Annotated Results
  # ============================================================
  
  message("==> Saving annotated results")
  start.time <- Sys.time()
  
  # Create output directory
  output_annotated <- file.path(output_folder, "TEs_annotated")
  
  tryCatch(
    dir.create(output_annotated, showWarnings = FALSE, recursive = TRUE),
    error = function(e) {
      stop(
        "Failed to create output directory '", output_annotated, "': ",
        e$message,
        call. = FALSE
      )
    }
  )
  
  # Save annotated results
  output.TE.res <- file.path(output_annotated, "DESeq2_TE_results_annotated.tsv")
  .save_deseq2(res.TEs, output.TE.res, "TE annotated results")
  
  res.TEs.down <- res.TEs %>% dplyr::filter(.data$log2FoldChange < -minlfc, .data$padj < maxpadj)
  output.TE.res.down <- file.path(output_annotated, "DESeq2_TE_results_down_annotated.tsv")
  .save_deseq2(res.TEs.down, output.TE.res.down, "TE down annotated results")

  res.TEs.up <- res.TEs %>% dplyr::filter(.data$log2FoldChange > minlfc, .data$padj < maxpadj)
  output.TE.res.up <- file.path(output_annotated, "DESeq2_TE_results_up_annotated.tsv")
  .save_deseq2(res.TEs.up, output.TE.res.up, "TE up annotated results")
  
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished saving results (", round(duration[[1]], 2), " seconds)")
   
  #============================================================
  # Step 7: Generate Region Distribution Plots
  #============================================================
  
  message("==> Generating region distribution plots")
  start.time <- Sys.time()
  
  # Update TE_results with annotated data
  TE_results_annotated <- list(
    res.TEs = res.TEs,
    TE.count = TE.count,
    gene.count = gene.count,
    metadata = metadata
  )
  
  tryCatch(
    {
      graphTEregion(
        TE_results_annotated,
        device = device,
        output_folder = output_annotated,
        width = width,
        height = height,
        minCounts = minCounts,
        plot.title = plot.title
      )
    },
    error = function(e) {
      warning(
        "Failed to generate region plots: ", e$message,
        "\nAnnotated results have been saved but plots were not created.",
        call. = FALSE
      )
    }
  )
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units = "secs")
  message("       -> Finished generating region distribution plots (", round(duration[[1]], 2), " seconds)")
  
  TE_results_annotated
  
}
