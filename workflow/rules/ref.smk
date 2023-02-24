rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v1.21.4/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    wrapper:
        "v1.21.4/bio/reference/ensembl-annotation"

rule chr_annotation:
    input:
        "resources/genome.gtf",
    output:
        temp("resources/genome_chr.gtf"),
    log:
        "logs/chr_annotation.log",
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        awk 'BEGIN{{ FS=OFS="\t" }} /^#/ {{ print }}; /^#/ {{ next }}; {{ gsub(/^/, "chr", $1); print }}' {input} > {output}
        """

rule chr_genome:
    input:
        "resources/genome.fasta",
    output:
        temp("resources/genome_chr.fasta"),
    log:
        "logs/chr_genome.log",
    conda:
        "../envs/htseq.yaml"
    shell:
        """
        awk '/^>/ {{ gsub(/^>/, ">chr"); print; next }} {{ print }}' {input} > {output}
        """

rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.21.4/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v1.21.4/bio/bwa/index"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "v1.21.4/bio/star/index"


rule RSEM_index:
    input:
        reference_genome="resources/genome.fasta",
    output:
        seq="index/reference.seq",
        grp="index/reference.grp",
        ti="index/reference.ti",
        tfa="index/reference.transcripts.fa",
        idxfa="index/reference.idx.fa",
        n2g="index/reference.n2g.idx.fa",
    params:
        extra="--gtf resources/genome.gtf",
    log:
        "logs/rsem/prepare-reference.log",
    wrapper:
        "v1.23.4/bio/rsem/prepare-reference"