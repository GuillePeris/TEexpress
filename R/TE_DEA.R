#' @import stringr
#' @title Compute differential expression analysis for transposable elements 
#'
#' @param metafile Full path of metadata file
#' @param folder Full path for folder of count files
#' @param output Path for output folder
#' @param useCtrlGenes Boolean for using only genes in estimating DESeq2 size factors (default: FALSE)
#' @param shrinklog2FC Boolean for applying log2FC shrinking in DESEq2 (default: FALSE)
#'
#' @returns Table of results
#' @export
#'
TE_DEA <- function(metafile, 
                   folder,
                   output=".",
                   useCtrlGenes=FALSE, 
                   shrinklog2FC=FALSE) {
  
  # Read metadata
  message("Reading metadata file.")
  metadata <- read_metadata(metafile)
  
  # Read counts into a count data frame
  message("Reading count files.")
  countData <- readTEcounts(metadata, folder)
  
  # DESeq2 analysis
  message("Calling DESeq2")
  dds <- call_deseq2(countData, metadata, useCtrlGenes)
  res <- results_deseq2(dds, shrinklog2FC)
  
  # Save normalized counts
  norm.counts <- norm_counts(dds)
  
  # Divide counts in genes and TEs and save
  gene.count <- norm.counts[!stringr::str_detect(rownames(norm.counts), ":"), ]
  output.gene<-paste0(output, "/",  "gene_normalizedCounts.csv")
  write.table(gene.count, file=output.gene, sep="\t")
  
  TE.count <- norm.counts[stringr::str_detect(rownames(norm.counts), ":"), ]
  output.TE<-paste0(output, "/",  "TE_normalizedCounts.csv")
  write.table(TE.count, file=output.TE, sep="\t")
  
  # Graphs
  
  
  # This is a procedure not returning anything
  invisible()
}
