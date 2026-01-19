#' Perform Differential Expression Analysis for Transposable Elements
#'
#' Main pipeline function for analyzing differential expression of transposable
#' elements (TEs) and genes using DESeq2. Reads count data, performs statistical
#' analysis, and generates visualization plots.
#'
#' @param metafile Character string. Full path to the metadata file. Should be
#'   a tab-separated file readable by \code{\link{read_metadata}}.
#' @param folder Character string. Full path to the directory containing count
#'   files listed in the metadata.
#' @param output Character string. Path for output directory. Will be created
#'   if it doesn't exist. Default is current directory (".").
#' @param maxpadj Numeric. Adjusted p-value threshold for significance.
#'   Features with padj < maxpadj are considered significant. Default is 0.05.
#' @param minlfc Numeric. Minimum absolute log2 fold change threshold for
#'   differential expression. Default is 1.
#' @param gtf.TE.file Character string. Full path to GTF file containing TE
#'   annotations. Must be the same GTF file used for upstream TE counting
#'   with TElocal from TEtranscript package.
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "svg", "eps", "png", "tiff", "jpeg". Can specify multiple formats
#'   as a vector (e.g., \code{c("jpeg", "png")}). Default is "png".
#' @param plot.title Character string. Title for volcano and MA plots. If empty
#'   string (default), generic titles will be used.
#' @param useCtrlGenes Logical. If TRUE, uses only genes (not TEs) for
#'   estimating DESeq2 size factors. This can improve normalization when TEs
#'   have very different expression patterns than genes. Default is FALSE.
#' @param shrinklog2FC Logical. If TRUE, applies apeglm log2 fold change
#'   shrinkage in DESeq2 for more accurate effect size estimates. Recommended
#'   for ranking and visualization. Default is FALSE.
#' @param saveNorm Logical. If TRUE, saves normalized count matrices for both
#'   genes and TEs. If FALSE, only saves DESeq2 results. Default is TRUE.
#'
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Reads sample metadata and count files
#'   \item Imports TE annotations from GTF file
#'   \item Adds genomic coordinates to TE features
#'   \item Filters non-standard chromosomes
#'   \item Performs DESeq2 differential expression analysis
#'   \item Extracts normalized counts
#'   \item Separates results for genes vs. TEs
#'   \item Saves results to separate directories
#'   \item Generates volcano and MA plots for both genes and TEs
#' }
#'
#' The function distinguishes TEs from genes by the presence of a colon (:) in
#' the feature name, which is the standard format for TE loci identifiers
#' (e.g., "L1Md_A:chr1:12345-12678:+").
#'
#' @section Output Structure:
#' The function creates the following directory structure:
#' \preformatted{
#' output/
#'   ├── genes_DEA/
#'   │   ├── DESeq2_gene_results.tsv
#'   │   ├── gene_normalizedCounts.tsv (if saveNorm = TRUE)
#'   │   ├── volcanoPlot.<format>
#'   │   └── maPlot.<format>
#'   └── TEs_DEA/
#'       ├── DESeq2_TE_results.tsv
#'       ├── TE_normalizedCounts.tsv (if saveNorm = TRUE)
#'       ├── volcanoPlot.<format>
#'       └── maPlot.<format>
#' }
#'
#' TE results include additional columns for genomic coordinates (seqnames,
#' start, end, strand, width) and TE class/family information.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{res.TEs}{Data frame of DESeq2 results for TEs, including genomic
#'     coordinates and TE annotations}
#'   \item{TE.count}{Data frame of normalized counts for TEs across samples}
#'   \item{gene.count}{Data frame of normalized counts for genes across samples}
#'   \item{metadata}{Original metadata data frame used for the analysis}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' results <- TE_DEA(
#'   metafile = "metadata.txt",
#'   folder = "counts",
#'   gtf.TE.file = "TEs.gtf"
#' )
#'
#' # Full analysis with custom parameters
#' results <- TE_DEA(
#'   metafile = "samples.txt",
#'   folder = "count_data",
#'   output = "results/DE_analysis",
#'   maxpadj = 0.01,
#'   minlfc = 1.5,
#'   gtf.TE.file = "annotations/TE_annotation.gtf",
#'   device = c("tiff", "png"),
#'   plot.title = "KO vs WT",
#'   useCtrlGenes = TRUE,
#'   shrinklog2FC = TRUE,
#'   saveNorm = TRUE
#' )
#'
#' # Access TE results
#' head(results$res.TEs)
#' dim(results$TE.count)
#'
#' # Find significantly upregulated TEs
#' sig_up_TEs <- results$res.TEs[
#'   results$res.TEs$padj < 0.05 &
#'   results$res.TEs$log2FoldChange > 1,
#' ]
#' }
#'
TE_DEA <- function(metafile,
                   folder,
                   output = ".",
                   maxpadj = 0.05,
                   minlfc = 1,
                   gtf.TE.file,
                   device = "png",
                   plot.title = "",
                   useCtrlGenes = FALSE,
                   shrinklog2FC = FALSE,
                   saveNorm = TRUE) {
  message("============================================") 
  message("  TE loci differential expression analysis  ")
  message("============================================") 
  
  # Input validation
  if (missing(metafile)) {
    stop("Argument 'metafile' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(metafile)) {
    stop("Metadata file '", metafile, "' not found.", call. = FALSE)
  }
  
  if (missing(folder)) {
    stop("Argument 'folder' is missing with no default.", call. = FALSE)
  }
  
  if (!dir.exists(folder)) {
    stop("Count data directory '", folder, "' not found.", call. = FALSE)
  }
  
  # Validate gtf.TE.file
  if (missing(gtf.TE.file)) {
    stop("Argument 'gtf.TE.file' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(gtf.TE.file)) {
    stop("TE GTF file '", gtf.TE.file, "' not found.", call. = FALSE)
  }
  
  for (argument in c(useCtrlGenes, shrinklog2FC, saveNorm)) {
    if (!is.logical(argument)) {
      stop("'", argument, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  # ============================================================
  # Step 1: Read Metadata and Count Files
  # ============================================================
  message("==> Reading metadata and count files.")
  start.time <- Sys.time()
  metadata <- tryCatch(
    read_metadata(metafile),
    error = function(e) {
      stop("Failed to read metadata: ", e$message, call. = FALSE)
    }
  )
  
  countData <- tryCatch(
    readTEcounts(metadata, folder),
    error = function(e) {
      stop("Failed to read count files: ", e$message, call. = FALSE)
    }
  )  
  
  if (ncol(countData) != nrow(metadata)) {
    stop(
      "Column count mismatch: count data has ", ncol(countData),
      " columns but metadata has ", nrow(metadata), " samples.",
      call. = FALSE
    )
  }
  
  end.time <- Sys.time()
  
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished reading metadata and count files (", 
          round(duration[[1]], 2), " seconds).")
  
  # ============================================================
  # Step 2: Load TE GTF Annotation
  # ============================================================
  
  message("==> Reading TE GTF file.")
  start.time <- Sys.time()
  
  gtf.TE <- tryCatch(
    importGTF(gtf.TE.file, format = "gtf"),
    error = function(e) {
      stop("Failed to import TE GTF file: ", e$message, call. = FALSE)
    }
  )  
  
  TE_annot.df <- as.data.frame(gtf.TE)
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished reading TE GTF file (", 
          round(duration[[1]], 2), " seconds).")  

  # ============================================================
  # Step 3: Add TE Coordinates and Filter
  # ============================================================
  
  message("==> Adding TE coordinates.")
  start.time <- Sys.time()
  
  # Separate TEs from genes based on presence of colon
  # TE format: "L1Md_A:chr1:12345-12678:+"
  # Gene format: "ENSMUSG00000000001"
  has_colon <- stringr::str_detect(rownames(countData), ":")
  
  countData.TEs <- countData[has_colon, ]
  countData.genes <- countData[!has_colon, ]
  
  # Add coordinates to TEs
  countData.TEs <- tryCatch(
    TE_coords(countData.TEs, TE_annot.df),
    error = function(e) {
      stop("Failed to add TE coordinates: ", e$message, call. = FALSE)
    }
  )
  
  # Remove TEs in non-standard chromosomes
  countData.TEs <- countData.TEs %>% 
     dplyr::filter(!is.na(.data$seqnames))
  
  # Recombine TEs and genes
  countData <- rbind(countData.TEs %>% 
                       dplyr::select(metadata$Sample), 
                     countData.genes)
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished adding TE coordinates (", 
          round(duration[[1]], 2), " seconds).") 
  
  # ============================================================
  # Step 4: Run DESeq2 Analysis
  # ============================================================  
  
  message("==> Running DESeq2 differential expression analysis")  
  if (useCtrlGenes) {
    message("    Using genes only for size factor estimation")
  }
  if (shrinklog2FC) {
    message("    Applying log2 fold change shrinkage")
  }
  start.time <- Sys.time()
  
  dds <- tryCatch(
    call_deseq2(countData, metadata, useCtrlGenes),
    error = function(e) {
      stop("DESeq2 analysis failed: ", e$message, call. = FALSE)
    }
  )
  
  res <- tryCatch(
    results_deseq2(dds, shrinklog2FC),
    error = function(e) {
      stop("Failed to extract DESeq2 results: ", e$message, call. = FALSE)
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished DESeq2 analysys (", 
          round(duration[[1]], 2), " seconds).") 
  
  # ============================================================
  # Step 5: Get Normalized Counts
  # ============================================================
  message("==> Extracting normalized counts")
  start.time <- Sys.time()
  
  norm.counts <- tryCatch(
    norm_counts(dds),
    error = function(e) {
      stop("Failed to get normalized counts: ", e$message, call. = FALSE)
    }
  )  
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished extracting normalized counts (", 
          round(duration[[1]], 2), " seconds).")
  
  # ============================================================
  # Step 6: Separate and Process Gene/TE Results  
  # ============================================================
  message("==> Separating and processing gene and TE results")  
  start.time <- Sys.time()
  
  # Separate gene and TE normalized counts 
  # Separate normalized counts
  has_colon_norm <- stringr::str_detect(rownames(norm.counts), ":")
  gene.count <- norm.counts[!has_colon_norm, ]
  TE.count <- norm.counts[has_colon_norm, ]
  
  # Separate DESeq2 results
  has_colon_res <- stringr::str_detect(rownames(res), ":")
  res.genes <- res[!has_colon_res, ]
  res.TEs <- res[has_colon_res, ]

  # Add TE metadata columns
  TE.count <- tryCatch(
    addTEColumns(TE.count, countData.TEs),
    error = function(e) {
      stop("Failed to add TE columns to normalized counts: ", e$message,
           call. = FALSE)
    }
  )
  
  res.TEs <- tryCatch(
    addTEColumns(res.TEs, countData.TEs),
    error = function(e) {
      stop("Failed to add TE columns to results: ", e$message, call. = FALSE)
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished processing gene and TE results (", 
          round(duration[[1]], 2), " seconds).")
  
  # ============================================================
  # Step 7: Save Results to Files
  # ============================================================
  
  message("==> Saving results to files")
  start.time <- Sys.time()
  
  output.genes <- paste0(output, "/genes_DEA")
  output.TEs <- paste0(output, "/TEs_DEA")

  tryCatch(
    {
      dir.create(output.genes, showWarnings = FALSE, recursive = TRUE)
      dir.create(output.TEs, showWarnings = FALSE, recursive = TRUE)
    },
    error = function(e) {
      stop("Failed to create output directories: ", e$message, call. = FALSE)
    }
  )
  
  # Save normalized counts if requested
  if (saveNorm) {
    output.genefile <- file.path(output.genes, "gene_normalizedCounts.tsv")
    tryCatch(
      {
        utils::write.table(
          gene.count,
          file = output.genefile,
          sep = "\t",
          quote = FALSE,
          row.names = TRUE,
          col.names = NA
        )
      },
      error = function(e) {
        warning("Failed to save gene normalized counts: ", e$message,
                call. = FALSE)
      }
    )
    
    output.TEfile <- file.path(output.TEs, "TE_normalizedCounts.tsv")
    tryCatch(
      {
        utils::write.table(
          TE.count,
          file = output.TEfile,
          sep = "\t",
          quote = FALSE,
          row.names = TRUE,
          col.names = NA
        )
      },
      error = function(e) {
        warning("Failed to save TE normalized counts: ", e$message,
                call. = FALSE)
      }
    )
  }
  
  # Save DESeq2 results
  output.gene.res.file <- file.path(output.genes, "DESeq2_gene_results.tsv")
  tryCatch(
    {
      utils::write.table(
        res.genes,
        file = output.gene.res.file,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
      )
    },
    error = function(e) {
      stop("Failed to save gene results: ", e$message, call. = FALSE)
    }
  )
  
  output.TE.res.file <- file.path(output.TEs, "DESeq2_TE_results.tsv")
  tryCatch(
    {
      utils::write.table(
        res.TEs,
        file = output.TE.res.file,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
      )
    },
    error = function(e) {
      stop("Failed to save TE results: ", e$message, call. = FALSE)
    }
  )
  
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished saving results (", 
          round(duration[[1]], 2), " seconds).")
  
  # ============================================================
  # Step 8: Generate Plots
  # ============================================================
  
  message("==> Generating visualization plots")
  start.time <- Sys.time()
  
  # Generate plots for genes
  tryCatch(
    {
      graphTE_DEA(
        res.genes, maxpadj, minlfc, device,
        output.genes,
        plot.title = if (nchar(plot.title) > 0) {
          paste0(plot.title, " - Genes")
        } else {
          "Genes"
        }
      )
    },
    error = function(e) {
      warning("Failed to generate gene plots: ", e$message, call. = FALSE)
    }
  )
  
  # Generate plots for TEs
  tryCatch(
    {
      graphTE_DEA(
        res.TEs, maxpadj, minlfc, device,
        output.TEs,
        plot.title = if (nchar(plot.title) > 0) {
          paste0(plot.title, " - TEs")
        } else {
          "TEs"
        }
      )
    },
    error = function(e) {
      warning("Failed to generate TE plots: ", e$message, call. = FALSE)
    }
  )
  end.time <- Sys.time()
  duration <- difftime(end.time, start.time, units="secs")
  message("       -> Finished generating plots (", 
          round(duration[[1]], 2), " seconds).")
  
  # Create return list
  TE_results <- list(
    res.TEs = res.TEs,
    TE.count = TE.count,
    gene.count = as.data.frame(gene.count),
    metadata = metadata
  ) 
  
  TE_results
}
