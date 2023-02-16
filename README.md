# Snakemake workflow: fastq2bakR

This is an extension of the rna-seq-star-deseq2 workflow described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2). It is designed to take as input fastq files from nucleotide recoding RNA-seq data (e.g., [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582#:~:text=TimeLapse%2Dseq%20is%20a%20single%2Dmolecule%20approach%20to%20monitor%20transcriptome,apparent%20using%20standard%20RNA%2Dseq.), [SLAM-seq](https://www.nature.com/articles/nmeth.4435), etc.). The output includes much of what is produced by the rna-seq-star-deseq2 workflow, but includes the necessary input for [bakR](https://github.com/simonlabcode/bakR), a tool I developed to perform [differential kinetic analysis](https://www.biorxiv.org/content/10.1101/2022.09.02.505697v1) with NR-seq data.

This pipeline is still in active development