
#' @import rtracklayer
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import ChIPseeker


TE_genomic_context <- function(res.TEs, gtf.genes,
                               TSSminus=-5000, TSSplus=5000,
                               downstream=10000) {
  # Setting ghost variables to NULL to pass check() 
  annotation <- NULL
  
  # Keep only canonical protein-coding genes
  gtf.genes <- gtf.genes[!is.na(gtf.genes$transcript_biotype) &
                             gtf.genes$transcript_biotype == "protein_coding", ]
  
  # Chipseeker options.
  options(ChIPseeker.downstreamDistance = downstream)
  options(ChIPseeker.ignore_1st_exon = TRUE)
  options(ChIPseeker.ignore_1st_intron = TRUE)
  
  TE.gr <- GenomicRanges::makeGRangesFromDataFrame(res.TEs)

  
  # Make TxDb for ChipSeeker
  gene.TxDb <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gtf.genes))

  # Annotate TEs to protein-coding genes with ChipSeeker
  anno <- ChIPseeker::annotatePeak(TE.gr, tssRegion=c(TSSminus, TSSplus), TxDb=gene.TxDb,
                       # Priorize coding regions
                       genomicAnnotationPriority=c("5UTR","3UTR",
                                                   "Exon","Promoter", "Intron",
                                                   "Downstream",
                                                   "Intergenic"),
                       ignoreDownstream = TRUE,
                       # overlap = "all" annotates with actual gene and not near TSS.
                       # Also, only consider promotors if overlapping TSS.
                       overlap = "all") 
  
  # annotatePeak may change internally UCSC notation to NCBI.
  anno@anno$geneChr <- ifelse(
    anno@anno$geneChr %in% c("X", "Y", "M", "MT"),
    paste0("chr", anno@anno$geneChr),
    ifelse(grepl("^[0-9]+$", anno@anno$geneChr),
           paste0("chr", anno@anno$geneChr),
           anno@anno$geneChr)
  )
  
  res.TEs <- cbind(res.TEs, mcols(anno@anno))
  res.TEs <- modifyAnnot(res.TEs)
  
  # Assign intergenic TEs to closest transcript
  transcript.gr <- gtf.genes[gtf.genes$type == "transcript", ]
  res.TEs.intergenic <- res.TEs %>% 
                           dplyr::filter(annotation == "Intergenic") %>% 
                           dplyr::select(seqnames, start, end, strand)
  TEs.intergenic.gr <- GenomicRanges::makeGRangesFromDataFrame(res.TEs.intergenic)
  intergenic_index <- GenomicRanges::nearest(TEs.intergenic.gr, transcript.gr)
  query <- transcript.gr[intergenic_index]
  my.columns <- c("geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneId")
  res.TEs[res.TEs$annotation == "Intergenic", my.columns] <- data.frame(as.character(seqnames(query)),
                                                              start(query),
                                                              end(query),
                                                              end(query)- start(query),
                                                              as.character(strand(query)),
                                                              query$gene_id)
  # Modify geneStrand notation
  res.TEs[!is.na(res.TEs$geneStrand) & res.TEs$geneStrand == '1', ]$geneStrand <-  '+'
  res.TEs[!is.na(res.TEs$geneStrand) & res.TEs$geneStrand == '2', ]$geneStrand <-  '-' 
  
  res.TEs
}

modifyAnnot <- function(df) {
  df$annotation[startsWith(df$annotation, "Intron")] <- "Intron"
  df$annotation[startsWith(df$annotation, "Exon")] <- "Exon"
  df$annotation[startsWith(df$annotation, "Distal")] <- "Intergenic"
  df$annotation[startsWith(df$annotation, "Downstream")] <- "Downstream"
  df$annotation[startsWith(df$annotation, "Promoter")] <- "Promoter"
  df
}
