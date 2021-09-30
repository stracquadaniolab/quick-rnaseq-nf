#!/usr/bin/env Rscript

'differential_expression.R

Usage:
  differential_expression.R <inputfile> <outputfile> [--control=<contrast_control>] [--case=<contrast_case>]

Options:
  --control=<contrast_control>                Condition to use as control [default: control].
  --case=<contrast_case>                      Condition to use as case [default: case].
  -h --help                                   Show this screen.
  --version                                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'quant_qc.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)

# performing differential expression analysis
results <- results(dse, contrast = c("condition", arguments$case, arguments$control))

# save results
write_csv(as_tibble(results), file=arguments$outputfile)

