#!/usr/bin/env Rscript

"quick-rnaseq-summarize-to-gene.R

Usage:
  quick-rnaseq-summarize-to-gene.R <samplesheet> [--inputdir=<inputdir>] [--output=<outfile>] [--counts-from-abundance=<counts_model>] [--db-organism=<db_org>]

Options:
  -i --inputdir=<inputdir>                    Directory containing salmon results [default: .].
  -o --output=<outfile>                       RDS output file [default: summarized-experiment.rds].
  -c --counts-from-abundance=<counts_model>   Generate estimated counts using abundance estimates either:
                                              no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
  -d --db-organism=<db_org>                   Organism database for gene annotation [default: org.Hs.eg.db]
  -h --help                                   Show this screen.
  --version                                   Show version.
" -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = "quick-rnaseq-summarize-to-gene.R")

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(tximeta))
suppressMessages(library(DESeq2))
suppressMessages(library(arguments$db_organism, character.only = TRUE))

# reading sample sheet
sample_sheet <- read_csv(arguments$samplesheet)

# prepare metadata for tximeta and deseq2
design <- tibble(
  names = sample_sheet$sample,
  files = file.path(arguments$inputdir, sample_sheet$sample, "quant.sf"),
  condition = factor(sample_sheet$condition)
)

# quantify and annotate transcripts
se <- tximeta(design, type = "salmon")

# summarize counts to gene level
gse <- summarizeToGene(se, countsFromAbundance = arguments$counts_from_abundance)

# adding gene symbols
gse <- addIds(gse, "SYMBOL", gene = TRUE)

# build deseq dataset
dds <- DESeqDataSet(gse, ~condition)

# fitting deseq model
dse <- DESeq(dds)

# save model
saveRDS(dse, file = arguments$output)