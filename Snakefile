# Snakefile

#configfile: "config.yaml"

rule all:
    input:
        "results/ML_rooted_tree.svg",
        "results/ML_final_tree.tiff" 

rule extract_clean:
    input:
        "data/seqdump.txt"
    output:
        fasta="results/cleaned.fasta",
        metadata="results/metadata.csv"
    conda:
        "environment.yml"
    shell:
        "python scripts/extract_and_clean.py"

rule align_mafft:
    input:
        "results/cleaned.fasta"
    output:
        "results/aligned.fasta"
    conda:
        "environment.yml"
    shell:
        "bash align/seq_align.bash {input} {output}"

rule iqtree_build:
    input:
        "results/aligned.fasta"
    output:
        "results/aligned.fasta.treefile"
    conda:
        "environment.yml"
    shell:
        "iqtree2 -s {input} -m GTR+G -bb 1000 -alrt 1000 -nt AUTO"

rule visualise_tree:
    input:
        tree="results/aligned.fasta.treefile",
        meta="results/metadata.csv"
    output:
        svg="results/ML_rooted_tree.svg"
    conda:
        "environment.yml"
    shell:
        "python scripts/visualise_tree.py"

rule final_R_plot:
    input:
        "results/aligned.fasta.treefile"
    output:
        "results/ML_final_tree.tiff"
    conda:
        "environment.yml"
    shell:
        "Rscript scripts/phylo_R.R"
