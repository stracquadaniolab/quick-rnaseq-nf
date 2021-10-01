#!/usr/bin/env Rscript

'quick-rnaseq-pca.R

Usage:
  quick-rnaseq-pca.R <inputfile> <outputfile>

Options:
  -h --help                                   Show this screen.
  --version                                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'quick-rnaseq-pca.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)

# code adapted from the deseq2 vignette
pca <- plotPCA(rlog(dse), intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca, "percentVar"))
plt <- ggplot(pca, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",pct_var[1],"% variance")) +
        ylab(paste0("PC2: ",pct_var[2],"% variance")) + 
        theme(panel.background = element_rect(fill = "white", colour = "#737373", size=1), legend.key=element_blank())

# saving to file with 4/3 ratio
ggsave(arguments$outputfile, plot = plt, width=5, height=5*(3/4))
