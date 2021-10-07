#!/usr/bin/env Rscript

"quick-rnaseq-cook.R

Usage:
  quick-rnaseq-cook.R <inputfile> <outputfile>

Options:
  -h --help                                   Show this screen.
  --version                                   Show version.
" -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = "quick-rnaseq-cook.R")

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

# reading deseq object
dse <- readRDS(arguments$inputfile)
cook_dist <- as_tibble(log10(assays(dse)[["cooks"]]), rownames = "gene_id") %>%
  pivot_longer(!gene_id)

# plotting cook distance by sample
plt <- ggplot(cook_dist, aes(x = name, y = value)) +
  geom_boxplot() +
  scale_y_continuous("log10 cook distance") +
  scale_x_discrete("") +
  theme(
    panel.background = element_rect(fill = "white", colour = "#737373", size = 1),
    legend.key = element_blank()
  )


# saving to file with 4/3 ratio
ggsave(arguments$outputfile, plot = plt, width = 6, height = 6 * (3 / 4))