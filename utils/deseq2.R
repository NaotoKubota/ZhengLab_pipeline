library(DESeq2)

args <- commandArgs(trailingOnly = T)
experiment_table <- args[1]
counts_table <- args[2]
output <- args[3]

print(paste("Experiment table:", experiment_table))
experiment <- read.table(

    experiment_table,
    sep = "\t",
    head = T,
    check.names = FALSE

)

print(paste("Counts table:", counts_table))
counts <- read.table(

    counts_table,
    sep = "\t",
    head = T,
    row.names = "PeakID",
    check.names = FALSE,
    stringsAsFactors = FALSE

)
counts_t <- t(counts)

print("Merge tables...")
merged <- merge(counts_t, experiment, by.x = 0, by.y = "bam")
print("Transpose...")
group <- data.frame(con = factor(merged$group))
counts <- merged[, -which (colnames(merged) %in% c("sample", "bam", "peak", "group"))]
rownames(counts) <- counts$Row.names
counts <- counts[, colnames(counts) != "Row.names"]
counts <- t(counts)
counts <- as.matrix(counts)

# At least six counts
print("At least six counts...")
counts <- counts[apply(counts, 1, sum)>6,]

# DEseq
print("DEseq2 analysis...")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = group, design = ~ con)
dds$con <- relevel(dds$con, ref = "Ref")
dds <- DESeq(dds)
res <- results(dds)
res$PeakID <- row.names(res)
res <- res[, c("PeakID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res <- res[order(res$padj), ]

# save
print(paste("Save results to:", output))
write.table(

    res,
    file = output,
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F

)
