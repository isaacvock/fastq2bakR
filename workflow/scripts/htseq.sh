#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
annotation=$4
output=$5
output2=$6
mutcnt=$7

    # Will create the ./results/htseq
    touch "$output2"


    samtools view -h -@ "$cpus" "$input" \
        | parallel \
            -L 2 \
            -j "$cpus" \
            --linebuffer \
            --joblog "$sample"_htseq_parallel.log \
            --roundrobin \
            --header '(@.*\n)*' \
            --pipe python $mutcnt \
                        -f sam \
                        --samout ./results/htseq/"$sample"_htseq.{#}_temp.sam \
                        -t transcript,exon,exon \
                        -i gene_id,gene_id,gene_id \
                        -m union,union,intersection-strict \
                        -c ./results/htseq/"$sample"_GF_htseq.{#}_temp.txt,./results/htseq/"$sample"_EF_htseq.{#}_temp.txt,./results/htseq/"$sample"_XF_htseq.{#}_temp.txt \
                        - \
                        "$annotation"

    # Error catching
    if [ $(cut -f 7 "$sample"_htseq_parallel.log | grep '1' | wc -l ) -eq 0 ]; then
        echo '* HTSeq counting finished for sample' $sample
    else
        echo '!!! HTSeq counting failed!'
        exit 1
    fi



# Combine outputs from individual jobs
    # .sam file
    samtools view -H "$input" \
        | cat - ./results/htseq/"$sample"_htseq.*_temp.sam \
        | samtools sort \
            -@ "$cpus" \
            -n \
            -o "$output" -

    echo "* HTSeq .sam files merged for sample $sample"

### Need to make this a separate rule!
    # gene counts
    parallel -j "$cpus" "awk -v 'OFS=\t' '
                                          {
                                                gene[\$1] += \$2
                                          }
                                          END {
                                                for (name in gene) {
                                                    print name, gene[name]
                                                }
                                          }' ./results/htseq/${sample}_{1}_htseq.*_temp.txt  \
                            | LC_COLLATE=C sort -k1,1 > ./results/htseq/{2}.${sample}_htseq.txt" ::: GF EF XF \
                                                                                 :::+ gene exon_w_overlaps mature


    echo "* HTSeq counts files merged for sample $sample"

    mv ${sample}_htseq_parallel.log ./results/htseq/
    rm ./results/htseq/${sample}*temp*