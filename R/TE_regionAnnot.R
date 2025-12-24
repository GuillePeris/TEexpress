#' @import GenomeInfoDb

TE_regionAnnot <- function(res.TEs, gtf.TE.file, gtf.genes.file, output_folder){

  # Load gtf's and check seqlevels 
  # BUSCAR FUNCIÓN PARA COMPROBAR GTFs EN TE_DENSITY
  gtf.TE <- formatSeqLevels(gtf.TE.file, format="gff")
  gtf.genes <- formatSeqLevels(gtf.genes.file, format="gff")
  TE_annot.df <- as.data.frame(gtf.TE)
  
  # Get coordinates from gtf.TE and add them to data frame
  res.TEs <- TE_coords(res.TEs, TE_annot.df)
  
  # Locate TEs relative to protein-coding genes
  res.TEs <- TE_genomic_context(res.TEs, gtf.genes)
  
  # Print annotated results
  output_folder <- paste0(output_folder, "/TEs_annotated")
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  output.TE.res<-paste0(output_folder, "/",  "DESeq2_TE_results_annotated.csv")
  write.table(res.TEs, file=output.TE.res, sep="\t")
  
  res.TEs
}

formatSeqLevels <- function(gtf.file, format = "gff") {
  gtf <- rtracklayer::import(gtf.file, format = format)
  gtf <- keepStandardChromosomes(gtf, pruning.mode="coarse")
  seqlevelsStyle(gtf) <- "UCSC"
  gtf
}
