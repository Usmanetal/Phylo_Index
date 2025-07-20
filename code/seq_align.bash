#!/bin/bash

# Alignment
mafft --auto results/cleaned.fasta > results/aligned.fasta

# Maximum‑likelihood tree

iqtree2 -s results/aligned.fasta -m GTR+G -bb 1000 -alrt 1000 -nt AUTO
# IQ‑TREE writes results/aligned.fasta.treefile  (plus .log, .iqtree, .contree …)


