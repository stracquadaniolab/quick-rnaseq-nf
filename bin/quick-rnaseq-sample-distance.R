#!/usr/bin/env Rscript

'quick-rnaseq-sample-distance.R

Usage:
  quick-rnaseq-sample-distance.R <inputfile> <outputfile> [--transform=<tf>]

Options:
  --transform=<tf>                            Transform function applied to counts [default: rlog]
  -h --help                                   Show this screen.
  --version                                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'quick-rnaseq-sample-distance.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)
dse_t <- NULL

# pick counts transformation
if (arguments$transform == "vst"){
  dse_t = vst(dse)
}else{
  dse_t = rlog(dse)
}

# compute the distance matrix
dist_matrix <- as.matrix(dist(t(assay(dse_t))))

# set pdf device with 4:3 ratio
pdf(arguments$outputfile, width=4, height=3)

# plot heatmap
heatmap(dist_matrix, labCol=FALSE, legend='col')

# close device
dev.off()


