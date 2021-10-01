#!/usr/bin/env Rscript

'qrna_go.R

Usage:
  qc_maplot.R <inputfile> <outputfile> [--fwer=<alpha>]

Options:
  -f --fwer=<alpha>           Family-wise error rate to filter differentially
                              expressed genes [default: 0.05].
  -d --db-organism=<db_org>   Organism database for gene annotation [default: org.Hs.eg.db]
  -g --gene-id=<gene_id>      Gene ID type [default: ensembl]
  -h --help                   Show this screen.
  --version                   Show version.
' -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = 'quant_qc.R')

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(GO.db))
suppressMessages(library(topGO))

# reading deseq object
res <- read_csv(arguments$inputfile)
res <- res[!is.na(res$pvalue),]

# gene list
gene_list <- res$pvalue
names(gene_list) <- res$gene_id

print(gene_list)

# GO data
go_data <- new("topGOdata",
                ontology = "BP",
                allGenes = gene_list,
                geneSelectionFun = function(.x) {.x <= (as.numeric(arguments$fwer)/length(gene_list))},
                annot = annFUN.org, mapping = arguments$db_organism, ID=arguments$gene_id)

# perform classical fisher test
fisher_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
fisher_results <- GenTable(go_data, classicFisher = fisher_test, topNodes = 1000)

# write results to file
write_csv(fisher_results, arguments$outputfile)
