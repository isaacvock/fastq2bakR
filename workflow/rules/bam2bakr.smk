
rule merge_bams:
    input:
        bamfiles=filter_bamfiles
    output:
        merged_bam="results/merge_bams/{sample}.merged.bam",
        sorted_bam="results/merge_bams/{sample}.sorted.bam"
    log:
        "logs/merge_bams/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "samtools merge -@ {threads} {output.merged_bam} {input.bamfiles} && samtools sort -n -@ {threads} {output.merged_bam} -o {output.sorted_bam} 1> {log} 2>&1"


rule sort_filter:
    input:
        "results/merge_bams/{sample}.sorted.bam",
    params:
        format = lambda wildcards: "PE" if get_format(wildcards) else "SE"
    output:
        "results/sf_reads/{sample}.s.sam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/sort_filter.sh {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1"

rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam",
        "resources/genome_chr.gtf",
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/htseq.sh {threads} {wildcards.sample} {input} {output} 1> {log} 2>&1"

rule normalize:
    input:
        expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
    output:
        "results/normalization/scale"
    params:
        spikename = config["spikename"]
    log:
        "logs/normalize/normalize.log"
    threads: 1
    conda:
        "../envs/normalize.yaml"
    shell:
        r"""
       	./workflow/scripts/normalize.R --dirs ./results/htseq/ --spikename {params.spikename}
	mv scale {output}
        """


rule call_snps:
    input:
        "resources/genome_chr.fasta",
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nsamps = nctl
    output:
        "results/snps/snp.txt",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/call_snps.sh {threads} {params.nsamps} {output} {input} 1> {log} 2>&1"

rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    params:
        format = lambda wildcards: "SE" if get_format(wildcards) else "PE",
        fragment_size = config["fragment_size"],
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = lambda wildcards: "F" if unique(get_strandedness(units)) == 'forward' else "R"
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/mut_call.sh {threads} {wildcards.sample} {input} {output} {params.fragment_size} {params.minqual} {params.mut_tracks} {params.format} {params.strand} 1> {log} 2>&1"

rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale",
        "resources/genome_chr.fasta",
    params:
        mut_tracks = config["mut_tracks"],
        wsl = config["WSL"],
        normalize = config["normalize"]
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    log:
        "logs/maketdf/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/tracks.yaml"
    shell:
        "./workflow/scripts/tracks.sh {threads} {wildcards.sample} {input} {params.mut_tracks} {params.wsl} {params.normalize} {output} 1> {log} 2>&1"

rule makecB:
    input:
        expand("results/counts/{sample}_counts.csv.gz", sample=SAMP_NAMES)
    params:
        mut_tracks = config["mut_tracks"],
        keepcols = config["keepcols"]
    output:
        "results/cB/cB.csv.gz"
    log:
        "logs/makecB/master.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "./workflow/scripts/master.sh {threads} {output} {params.keepcols} {params.mut_tracks} 1> {log} 2>&1"
