library(DESeq2)
setwd("C:/Users/danie/OneDrive/Documentos/GitHub/hungerGamesModel/files/strainSummaries/ri/genes")

library("tximport")
library("readr")
library("tximportData")


samples <- read.table("ri_samples_t04vst12.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('ri_counts_t04vst12.txt',sep="\t",row.names="patricID"))

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

samples <- read.table("ri_samples_t04vst48.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('ri_counts_t04vst48.txt',sep="\t",row.names="patricID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t04vst48_deseq.txt", sep='\t')

#############################################################################

samples <- read.table("ri_samples_t12vst48.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('ri_counts_t12vst48.txt',sep="\t",row.names="patricID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t12vst48_deseq.txt", sep='\t')


#############################################################################

samples <- read.table("btri_ri_samples_t12vst36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('btri_ri_counts_t12vst36.txt',sep="\t",row.names="patricID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="btri_ri_t12vst36_deseq.txt", sep='\t')


#############################################################################

samples <- read.table("btri_ri_samples_t12t48vst12t24t36.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('btri_ri_counts_t12t48vst12t24t36.txt',sep="\t",row.names="patricID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="btri__ri_t12t48vst12t24t36_deseq.txt", sep='\t')



#############################################################################

samples <- read.table("ri_samples_t12vsbtri_ri_t12.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('ri_counts_t12vsbtri_ri_t12.txt',sep="\t",row.names="patricID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.0)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="ri_t12vsbtri_ri_t12_deseq.txt", sep='\t')

