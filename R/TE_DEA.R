#' @import stringr
#' @importFrom utils write.table
#' @title Compute differential expression analysis for transposable elements 
#'
#' @param metafile Full path of metadata file
#' @param folder Full path for folder of count files
#' @param output Path for output folder (default = ".")
#' @param maxpadj P-adjusted value for significant features (default = 0.05)
#' @param minlfc Value for dysregulated features (default = 1.0)
#' @param device Format for graphs ("pdf", "svg", "eps", "png", "tiff", "jpeg"). 
#'               A vector for several formats can be uses: c("svg", "png") (Default: "png)
#' @param useCtrlGenes Boolean for using only genes in estimating DESeq2 size factors (default: FALSE)
#' @param shrinklog2FC Boolean for applying log2FC shrinking in DESEq2 (default: FALSE)
#'
#' @returns Table of results
#' @export
#'
TE_DEA <- function(metafile, 
                   folder,
                   output=".",
                   maxpadj = 0.05,
                   minlfc = 1,
                   device = "png",
                   useCtrlGenes=FALSE, 
                   shrinklog2FC=FALSE) {
  
  # Read metadata
  message("==> Reading metadata file.")
  metadata <- read_metadata(metafile)
  
  # Read counts into a count data frame
  message("==> Reading count files.")
  countData <- readTEcounts(metadata, folder)
  
  # DESeq2 analysis
  message("==> Calling DESeq2")
  dds <- call_deseq2(countData, metadata, useCtrlGenes)
  res <- results_deseq2(dds, shrinklog2FC)
  
  # Save normalized counts
  message("==> Saving normalized counts")
  norm.counts <- norm_counts(dds)

  # Create separate folders for genes and TEs
  output.genes <- paste0(output, "/genes")
  dir.create(output.genes, showWarnings = FALSE, recursive = TRUE)
  output.TEs <- paste0(output, "/TEs")
  dir.create(output.TEs, showWarnings = FALSE, recursive = TRUE)
    
  # Divide counts in genes and TEs and save
  gene.count <- norm.counts[!stringr::str_detect(rownames(norm.counts), ":"), ]
  output.gene<-paste0(output.genes, "/",  "gene_normalizedCounts.csv")
  write.table(gene.count, file=output.gene, sep="\t")
  
  TE.count <- norm.counts[stringr::str_detect(rownames(norm.counts), ":"), ]
  output.TE<-paste0(output.TEs, "/",  "TE_normalizedCounts.csv")
  write.table(TE.count, file=output.TE, sep="\t")
  
  # Graphs
  message("==> Creating graphs")
  
  res.genes <- res[!stringr::str_detect(rownames(res), ":"), ]
  graphTools(res.genes, maxpadj, minlfc, device, 
             output.genes, plot.title = "Prueba") 
  res.TEs <- res[stringr::str_detect(rownames(res), ":"), ]
  graphTools(res.TEs, maxpadj, minlfc, device, 
             output.TEs, plot.title = "Prueba") 
  
  # This is a procedure not returning anything
  invisible()
}
