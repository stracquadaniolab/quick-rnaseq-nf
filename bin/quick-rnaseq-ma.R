#!/usr/bin/env Rscript

'quick-rnaseq-ma.R

Usage:
  quick-rnaseq-ma.R <inputfile> <outputfile> [--control=<contrast_control>] [--case=<contrast_case>] [--transform=<tf>]

Options:
  --control=<contrast_control>                Condition to use as control [default: control].
  --case=<contrast_case>                      Condition to use as case [default: case].
  -h --help                                   Show this screen.
  --version                                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'quick-rnaseq-ma.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)
res <- lfcShrink(dse, contrast = c("condition", arguments$case, arguments$control), type="normal")
ma <- plotMA(res, returnData=TRUE)

# drawing a maplot using ggplot2
lfc <- ceiling(max(abs(ma$lfc), na.rm = T))
plt <- ggplot(ma, aes(x=mean, y=lfc, color=isDE)) + 
        geom_point() + 
        geom_hline(yintercept = 0, linetype = "dashed", colour = "#737373") + 
        scale_y_continuous("log fold change", limits = c(-lfc,lfc)) + 
        scale_x_continuous("mean of normalized counts") +  
        scale_colour_manual(values =c('#CCCCCC', '#08519c')) + 
        theme(panel.background = element_rect(fill = "white", colour = "#737373", size=1)) + guides(color = "none")
ggsave(arguments$outputfile, plot = plt, width=5, height=5*(3/4))
