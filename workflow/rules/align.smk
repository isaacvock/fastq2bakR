rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
    output:
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
        aln_tx="results/star/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        idx=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    conda:
        "../envs/star.yaml"
    threads: 24
    script: 
        "../scripts/star-align.py"


rule quantify:
    input:
        bam="results/star/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        reference=multiext("index/reference", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
    output:
        genes_results="results/rsem/{sample}-{unit}.genes.results",
        isoforms_results="results/rsem/{sample}-{unit}.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired-end=True,
        # additional optional parameters to pass to rsem, for example,
        extra="--seed 42",
    log:
        "logs/rsem/calculate_expression/{sample}-{unit}.log",
    conda:
        "../envs/rsem.yaml"
    threads: 2
    script:
        "../scripts/rsem-calc.py"