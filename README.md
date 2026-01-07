# TEexpress

Tools to analyze transposable elements expression, annotate them in genomic regions and classify them as expressed from their own promoter or from a runthrough transcription.

## Installation

You can install the development version of TEexpress from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("GuillePeris/TEexpress", 
                        upgrade = "never",
                        dependencies = TRUE)
```

## Requirements

This package has been designed explicitly for *TElocal* software from *TEtranscripts* package.

A tab-separated text file must be provided with this structure:
``` r
File Sample	Group	Condition
TElocal_PA1_WT1_chr22.cntTable	WT1	WT	Control
TElocal_PA1_WT2_chr22.cntTable	WT2	WT	Control
TElocal_PA1_WT3_chr22.cntTable	WT3	WT	Control
TElocal_PA1_WT4_chr22.cntTable	WT4	WT	Control
TElocal_PA1_KO1_chr22.cntTable	KO1	KO	Treat
TElocal_PA1_KO2_chr22.cntTable	KO2	KO	Treat
TElocal_PA1_KO3_chr22.cntTable	KO3	KO	Treat
TElocal_PA1_KO4_chr22.cntTable	KO4	KO	Treat
```
It is important that in the fourth column appears the tags "Control" and "Treat".

## TE differential expression analysis
First we need to perform differential expression analysis for TE counts. For that
purpose we have **TE_DEA** function.

```r
library(TEexpress)
metafile <- system.file("extdata", "data.csv", package = "TEexpress") # Full path to table file
folder <- system.file("extdata", package = "TEexpress")  # Full path to count file folder
output <- "results" # Full path to new results folder
maxpadj <- 0.05     # P-adjusted value for significant features (default = 0.05)
minlfc <- 1         # Value for dysregulated features (default = 1.0)
device <- c("png", "svg") # Format for graphs
gtf.TE.file <- ""
plot.title <- "Test PA1 DGCR8-KO versus WT" # Title for graphs
gtf.TE.file <- system.file("extdata", "Homo_sapiens.GRCh38.115_chr1.gtf", package = "TEexpress") # Full path to table file
TE_results <- TE_DEA(metafile, folder, output, maxpadj, minlfc, gtf.TE.file, device, plot.title)

```

## TE annotation to genomic regions
```r
gtf.genes.file <- system.file("extdata", "Homo_sapiens.GRCh38.115_chr1.gtf", package = "TEexpress") # Full path to table file
TE_results <- annotate_TE_regions(TE_results, gtf.genes.file, output, device, plot.title)
```

## TE classification as self-expressed or gene-dependent

To be done

