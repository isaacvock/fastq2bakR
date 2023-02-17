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
    threads: 24
    conda:
        "../envs/star.yaml"
    script: 
        "../scripts/star-align.py"