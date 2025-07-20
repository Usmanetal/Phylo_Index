#!/usr/bin/env python
"""
Roots the IQ‑TREE ML tree, renames tips with rich labels,
saves a circular SVG ready for publication.
"""

from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import pandas as pd
from pathlib import Path
from ete3 import Tree

TREEPATH = Path("results/aligned.fasta.treefile")
METAPATH = Path("results/metadata.csv")
OUTSVG   = Path("results/ML_rooted_tree.svg")

# --- load tree -------------------------------------------------
t = Tree(str(TREEPATH), format=1)   # Newick with branch lengths

# midpoint‑root (change to .set_outgroup() if you have a reference)
# t.set_outgroup(t.get_midpoint_outgroup())
t.set_outgroup(t & "HQ843607.1")



# --- load metadata --------------------------------------------
meta = pd.read_csv(METAPATH)
label_map = dict(zip(meta["accession"], meta["label"]))

# rename each leaf
for leaf in t:
    if leaf.name in label_map:
        leaf.name = label_map[leaf.name]

# --- style -----------------------------------------------------
for n in t.traverse():
    nstyle = NodeStyle()
    nstyle["size"] = 4
    nstyle["fgcolor"] = "black"
    n.set_style(nstyle)

ts = TreeStyle()
ts.mode = "r"               # circular
ts.show_leaf_name = True
ts.scale = 50
ts.title.add_face(TextFace("Rooted ML Tree – Accession | Country | Isolate | Year",
                           fsize=10, bold=True), column=0)

# export
OUTSVG.parent.mkdir(exist_ok=True, parents=True)
t.render(str(OUTSVG), w=900, tree_style=ts)
print(f"✅ SVG saved → {OUTSVG}")


# import pandas as pd
# import toytree
# import toyplot.svg
# import toyplot.pdf
# import toyplot.png

# # Load tree
# tree = toytree.tree("aligned.fasta.treefile")
# tree = tree.root("FN599684.1")  # Root on known tip

# # Load metadata
# meta = pd.read_csv("results/metadata.csv")


# print(tree.get_tip_labels()[:5])  # example tips
# print(meta['label'].head())       # sample metadata labels

# label_map = meta.set_index("accession")["label"].to_dict()



# # Map labels to colors using country
# tip_colors = meta.set_index("label")["country"].map(color_map).to_dict()


# # label_map is a dict: old_name -> new_label

# for tip in tree.tips():
#     if tip.name in label_map:
#         tip.name = label_map[tip.name]

import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo
import pandas as pd
from pathlib import Path

# --- Paths ---
TREEPATH = Path("results/aligned.fasta.treefile")
METAPATH = Path("results/metadata.csv")
OUTPDF   = Path("results/ML_rooted_tree.pdf")

# --- Load tree ---
tree = Phylo.read(str(TREEPATH), "newick")

# Remove inner node labels (like bootstrap values)
for clade in tree.get_nonterminals():
    clade.confidence = None
    clade.name = None

# --- Rename tip labels using metadata ---
meta = pd.read_csv(METAPATH)
label_map = dict(zip(meta["accession"], meta["label"]))

for leaf in tree.get_terminals():
    if leaf.name in label_map:
        leaf.name = label_map[leaf.name]

# --- Plot settings ---
fig = plt.figure(figsize=(10, 60), dpi=300)  # Taller figure for long trees
matplotlib.rc('font', size=6)                # Leaf & node labels
matplotlib.rc('xtick', labelsize=10)         # X ticks
matplotlib.rc('ytick', labelsize=10)         # Y ticks

# Create axes
ax = fig.add_subplot(1, 1, 1)

# Draw tree (no internal node labels)
Phylo.draw(tree, axes=ax, do_show=False, branch_labels=lambda c: None)

# Optional: Add a title
ax.set_title("Rooted ML Tree – Accession | Country | Isolate | Year", fontsize=12, weight='bold')

# Save the figure
OUTPDF.parent.mkdir(exist_ok=True, parents=True)
fig.savefig(OUTPDF, bbox_inches="tight")
print(f"✅ Tree plot saved → {OUTPDF}")
