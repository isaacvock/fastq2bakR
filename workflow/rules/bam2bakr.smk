import glob

import pandas as pd

unit_s4U = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name", "s4U"], drop=False)
    .sort_index()
)

names = list(unit_s4U.index.values)


SAMP_NAMES = list(zip(*names))[0]
UNIT_NAMES = list(zip(*names))[1]
S4U = list(zip(*names))[2]

CTL_NAMES = list([SAMP_NAMES[i] for i in range(len(SAMP_NAMES)) if S4U[i] == 'no'])

nctl = len(CTL_NAMES)

inputfiles = pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})

# Define a function to filter the input BAM files by sample name
def filter_bamfiles(wildcards):
    return "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, unit=inputfiles[inputfiles['sample_name'] == wildcards.sample]['unit_name'].values[0])

rule merge_bams:
    input:
        bamfiles=filter_bamfiles
    output:
        merged_bam="results/merge_bams/{sample}.merged.bam"
    log:
        "logs/merge_bams/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/merge.yaml"
    shell:
        "samtools merge -@ {threads} {output.merged_bam} {input.bamfiles}"


rule sort_filter:
    input:
        "results/merge_bams/{sample}.merged.bam",
    params:
        format="SE" if not is_paired_end(wildcards.sample) else "PE"
    output:
        "results/sf_reads/{sample}.s.sam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/sort.yaml"
    shell:
        "workflow/scripts/sort_filter.sh {threads} {wildcards.sample} {input} {output} {params.format}"

rule htseq_cnt:
    input:
        "results/sf_reads/{sample}.s.sam",
        "resources/genome.gtf",
    output:
        "results/htseq/{sample}_tl.bam",
        temp("results/htseq/{sample}_check.txt")
    log:
        "logs/htseq_cnt/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/htseq.yaml"
    shell:
        "workflow/scripts/htseq.sh {threads} {wildcards.sample} {input} {output}"

rule normalize:
    input:
        expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
    output:
        "results/normalization/scale"
    log:
        "logs/normalize/normalize.log"
    threads: 1
    conda:
        "../envs/normalize.yaml"
    shell:
        r"""
       	./workflow/scripts/normalize.R --dirs ./results/htseq/ --spikename {config[spikename]}
	mv scale {output}
        """


rule call_snps:
    input:
        "resources/genome.fasta",
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
        "../envs/snps.yaml"
    shell:
        "workflow/scripts/call_snps.sh {threads} {params.nsamps} {output} {input}"

rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    params:
        format="SE" if not is_paired_end(wildcards.sample) else "PE"
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: workflow.cores
    conda:
        "../envs/cnt_muts.yaml"
    shell:
        "workflow/scripts/mut_call.sh {threads} {wildcards.sample} {input} {output} {config[fragment_size]} {config[minqual]} {config[mut_tracks]} {params.format}"

rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale",
        "resources/genome.fasta",
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    threads: workflow.cores
    conda:
        "../envs/tracks.yaml"
    shell:
        "workflow/scripts/tracks.sh {threads} {wildcards.sample} {input} {config[mut_tracks]} {config[WSL]} {config[normalize]} {output}"

rule makecB:
    input:
        expand("results/counts/{sample}_counts.csv.gz", sample=config["samples"])
    output:
        "results/cB/cB.csv.gz"
    log:
        "logs/makecB/master.log"
    threads: workflow.cores
    conda:
        "../envs/cB.yaml"
    shell:
        "workflow/scripts/master.sh {threads} {output} {config[keepcols]} {config[mut_tracks]}"
