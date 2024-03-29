log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial
# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(data.frame("gene"=rownames(res),res), file=snakemake@output[["table"]], row.names=FALSE, sep='\t')
