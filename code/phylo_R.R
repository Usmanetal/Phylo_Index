
library(Biostrings)

library(DECIPHER)
library(ape)
library(phangorn)
library(ggtree)
library(treeio)
library(tidyverse)
library(readr)
meta <- read_csv("results/metadata.csv")
# 1. Read the tree correctly (creates a phylo/treedata object)
tr <- read.tree("results/aligned.fasta.treefile")
# Keep only 'accession' and 'label' (rename for plot clarity)
meta_clean <- as.data.frame(meta) %>%
  select(accession, label) %>%
  rename(tip_label = label)

# Build the plot
p <- ggtree(tr, layout = "rectangular") %<+% meta_clean +
     geom_tiplab(aes(label = tip_label), size = 2.5) +
     ggtitle("Rooted ML Tree with Custom Tip Labels") +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save the plot
ggsave("pdf/tree_with_labels.pdf", plot = p, width = 10, height = 20)
