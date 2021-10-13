#!/usr/bin/env Rscript

"quick-rnaseq-dge.R

Usage:
  quick-rnaseq-dge.R <inputfile> <outputfile> [--control=<contrast_control>] [--case=<contrast_case>] [--log-foldchange=<lfc>][--fdr=<alpha>]

Options:
  --control=<contrast_control>  Condition to use as control [default: control].
  --case=<contrast_case>        Condition to use as case [default: case].
  -l --log-foldchange=<lfc>     Log2 fold-change threshold [default: 0].
  -f --fdr=<alpha>              False Discovery Rate threshold for independent filtering [default: 0.01].
  -h --help                     Show this screen.
  --version                     Show version.
" -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = "quick-rnaseq-dge.R")

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)

# performing differential expression analysis
results <- results(dse,
  contrast = c("condition", arguments$case, arguments$control),
  saveCols = c("SYMBOL", "gene_id", "ENTREZID"),
  independentFiltering = TRUE,
  alpha = as.numeric(arguments$fdr),
  lfcThreshold = as.numeric(arguments$log_foldchange)
)

results_df <- as_tibble(results) %>%
  relocate(SYMBOL) %>%
  relocate(ENTREZID) %>%
  relocate(gene_id) %>%
  arrange(pvalue)

# save results
write_csv(as_tibble(results_df), file = arguments$outputfile)