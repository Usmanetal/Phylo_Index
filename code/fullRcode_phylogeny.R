# install.packages(c("ape", "phangorn", "BiocManager"))
# BiocManager::install(c("treeio", "Biostrings", "DECIPHER", "ggtree"))

library(ape)
library(phangorn)
library(treeio)
library(Biostrings)
library(DECIPHER)
library(ggtree)
library(tidyverse)

# Load the sequences
seqs <- readDNAStringSet("seqdump.txt", format = "fasta")


# Remove duplicates if any
seqs <- unique(seqs)

# Align using DECIPHER
aligned <- AlignSeqs(seqs, anchor = NA)
writeXStringSet(aligned, filepath = "aligned.fasta")

# Convert alignment to phangorn format
alignment <- read.phyDat("aligned.fasta", format = "fasta")

# Distance tree (initial guess)
dm <- dist.ml(alignment)
tree_init <- NJ(dm)
 
# ML optimization
fit <- pml(tree_init, data = alignment)
fit_opt <- optim.pml(fit, model = "GTR", optGamma = TRUE, rearrangement = "stochastic")

# Original tip labels
old_labels <- fit_opt$tree$tip.label

# Extract components: accession, country (if available), isolate
# new_labels <- sapply(old_labels, function(label) {
#   # Extract GenBank accession (first word)
#   accession <- sub(" .*", "", label)
  
#   # Try to extract isolate info (anything after "isolate" or "isolate ")
#   isolate <- sub(".*isolate[ ]*", "", label)
#   isolate <- sub("[^A-Za-z0-9_\\-].*", "", isolate)  # strip trailing junk

#   # Optionally extract country if it's embedded — this part is heuristic
#   country <- NA
#   if (grepl(" [A-Z][A-Za-z]+", label)) {
#     words <- unlist(strsplit(label, " "))
#     # Try matching country names heuristically (you can refine this)
#     known_countries <- c("USA", "Ghana", "Cameroon", "Guinea-Bissau", "Gabon", "Democratic Republic of the Congo", 
#     "Angola","Spain","Liberia", "Pakista", "Cyprus", "Sweden", "Guinea-Bissau",
#      "Togo", "Kenya", "Nigeria", "Uganda", "India", "Brazil", "China", "South Korea", "Russia")
#     match <- intersect(words, known_countries)
#     if (length(match) > 0) {
#       country <- match[1]
#     }
#   }

#   # Combine clean label
#   paste0(accession,
#          if (!is.na(country)) paste0("|", country) else "",
#          if (nchar(isolate) > 0) paste0("|", isolate) else "")
# })


new_labels <- sapply(old_labels, function(label) {
  # Extract GenBank accession number
  accession <- sub(" .*", "", label)

  # Extract isolate or clone name (allow letters, numbers, dash, underscore, dot)
  isolate <- NA
  if (grepl("isolate", label, ignore.case = TRUE)) {
    iso <- sub(".*isolate[ ]*", "", label, ignore.case = TRUE)
    isolate <- sub("[^A-Za-z0-9_.-].*", "", iso)
  } else if (grepl("clone", label, ignore.case = TRUE)) {
    clo <- sub(".*clone[ ]*", "", label, ignore.case = TRUE)
    isolate <- sub("[^A-Za-z0-9_.-].*", "", clo)
  }

  # Extract country (one or two capitalized words after 'from')
  country <- NA
  if (grepl("from [A-Z][a-z]+( [A-Z][a-z]+)?", label)) {
    country <- sub(".*from ([A-Z][a-z]+(?: [A-Z][a-z]+)?).*", "\\1", label)
  }

  # Combine into final cleaned label
  paste0(accession,
         if (!is.na(country)) paste0("|", country) else "",
         if (!is.na(isolate)) paste0("|", isolate) else "")
})



# Assign the new labels to the tree
fit_opt$tree$tip.label <- new_labels


# Root the tree (choose a known accession manually if applicable)
# rooted_tree <- root(fit_opt$tree, outgroup = "FN599698.1|255HALD", resolve.root = TRUE)
rooted_tree <- root(fit_opt$tree, outgroup = "HQ843646.1|Nigeria|07NG.SN131", resolve.root = TRUE)

# Save rooted tree
write.tree(rooted_tree, file = "rooted_tree.nwk")

# Optional metadata (if you want to add host/country later)
meta <- data.frame(label = rooted_tree$tip.label)

# Plot
p <- ggtree(rooted_tree, layout = "circular") +
  geom_tiplab(size = 2, aes(label = label), hjust = -0.1) +
  geom_tippoint(color = "steelblue", size = 2) +
  ggtitle("Rooted ML Phylogenetic Tree of HIV-1 pol Region") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save
ggsave("HIV1_phylo_tree.tiff", plot = p, width = 10, height = 10, dpi = 300)

tip_labels <- fit_opt$tree$tip.label
"FN599684" %in% tip_labels



tip_metadata <- data.frame(
  original_label = old_labels,
  cleaned_label = new_labels,
  accession = sub(" .*", "", old_labels),
  isolate = sapply(old_labels, function(x) {
    iso <- sub(".*isolate[ ]*", "", x, ignore.case = TRUE)
    sub("[^A-Za-z0-9_\\-].*", "", iso)
  }),
  country = sapply(old_labels, function(x) {
    if (grepl("from [A-Z][a-z]+", x)) {
      sub(".*from ([A-Z][a-z]+).*", "\\1", x)
    } else {
      NA
    }
  }),
  stringsAsFactors = FALSE
)

write.csv(tip_metadata, "tip_metadata.csv", row.names = FALSE)

rooted_tree <- ape::drop.tip(rooted_tree, tip = "KX398187.1|Cameroon|YBF282")
# Set terminal (tip) branch lengths to 0
rooted_tree$edge.length[rooted_tree$edge[,2] <= length(rooted_tree$tip.label)] <- 0


# Step 4: Plot tree with gheatmap
p1 <- ggtree(rooted_tree) +
  geom_tiplab(size = 2, align = TRUE)
p1
rooted_tree$edge.length <- rooted_tree$edge.length / 10  # or another factor


p2 <- ggtree(rooted_tree) +
  geom_treescale() +  # no x/y needed
  geom_tiplab(size = 2)

p2

ggsave("my_tree_plot.pdf", p2, width = 8, height = 10)
p2 <- ggtree(rooted_tree, layout = "rectangular") +
  geom_tiplab(size = 2, offset = 0.01) +  # smaller offset
  geom_treescale()
p2
p3 <- ggtree(rooted_tree) %<+% tip_metadata +
  geom_tiplab(aes(color = country), size = 2)

ggsave("my_tree_plot.png", plot = p2, width = 8, height = 10, dpi = 300)

library(ape)
aligned <- read.dna("aligned.fasta", format = "fasta")
dist_mat <- dist.dna(aligned)
heatmap(as.matrix(dist_mat), col = viridis::viridis(100))


# Load sequences
seqs <- readDNAStringSet("seqdump.txt")
names(seqs) <- make.names(names(seqs), unique = TRUE)

# Detect known DRMs in the header
known_drms <- c("M184V", "K65R", "Y181C", "K103N", "D67N", "T215Y")  # Add more as needed

extract_drm_status <- function(desc, drm_list) {
  any(sapply(drm_list, function(drm) grepl(drm, desc)))
}

metadata <- data.frame(
  label = names(seqs),
  description = names(seqs),
  drug_resistant = sapply(names(seqs), function(x) extract_drm_status(x, known_drms))
)

aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, "aligned.fasta")

# ML Tree
phydat <- read.phyDat("aligned.fasta", format = "fasta")
dm <- dist.ml(phydat)
tree_init <- NJ(dm)
fit <- pml(tree_init, data = phydat)
fit_opt <- optim.pml(fit, model = "GTR", optGamma = TRUE, rearrangement = "stochastic")
rooted_tree <- midpoint(fit_opt$tree)

p <- ggtree(rooted_tree) %<+% metadata +
  geom_tippoint(aes(color = drug_resistant), size = 2) +
  geom_tiplab(aes(label = label), size = 2, align = TRUE) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  ggtitle("Phylogenetic Tree Highlighting Drug-Resistant Sequences") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("tree_DRMs_colored.tiff", plot = p, width = 12, height = 14, dpi = 600)

