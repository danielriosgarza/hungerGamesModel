library(DESeq2)
setwd("C:/Users/danie/OneDrive/Documentos/GitHub/hungerGamesModel/files/strainSummaries/bt/genes")

library("tximport")
library("readr")
library("tximportData")


samples <- read.table("bt_samples_t04vst12.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t04vst12.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t04vst12_deseq.txt", sep='\t')

#############################################################################

samples <- read.table("bt_samples_t04vst36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t04vst36.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t04vst36_deseq.txt", sep='\t')


#############################################################################

samples <- read.table("bt_samples_t12vst36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t12vst36.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t12vst36_deseq.txt", sep='\t')


#############################################################################

samples <- read.table("bt_samples_t12vst36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t12vst36.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t12vst36_deseq.txt", sep='\t')


#############################################################################

samples <- read.table("bt_samples_t36vs_btri_bt_t36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t36vs_btri_bt_t36.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t36vs_btri_bt_t36_deseq.txt", sep='\t')


samples <- read.table("bt_samples_t12t36vs_btri_bt_t12t24t36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bt_counts_t12t36vs_btri_bt_t12t24t36.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t12t36vs_btri_bt_t12t24t36_deseq.txt", sep='\t')
