#!/usr/bin/env Rscript

'summarize_to_gene.R

Usage:
  summarize_to_gene.R <samplesheet> [--inputdir=<inputdir>] [--output=<outfile>] [--counts-from-abundance=<counts_model>]

Options:
  -i --inputdir=<inputdir>                    Directory containing salmon results [default: .].
  -o --output=<outfile>                       RDS output file [default: summarized-experiment.rds].
  -c --counts-from-abundance=<counts_model>   Generate estimated counts using abundance estimates either:
                                              no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
  -h --help                                   Show this screen.
  --version                                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'summarize_to_gene.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(tximeta))
suppressMessages(library(DESeq2))

# reading sample sheet
sample_sheet <- read_csv(arguments$samplesheet)

# prepare metadata for tximeta and deseq2
design <- tibble(names=sample_sheet$sample, files=file.path(arguments$inputdir,paste(sample_sheet$sample, "_quant", sep=''), "quant.sf"), condition = sample_sheet$condition)

# quantify and annotate transcripts
se <- tximeta(design)

# summarize counts to gene level
gse <- summarizeToGene(se, countsFromAbundance = arguments$counts_from_abundance)

# build deseq dataset
dds <- DESeqDataSet(gse, ~ condition)

# fitting deseq model
dse <- DESeq(dds)

# save model
saveRDS(dse,file = arguments$output)

