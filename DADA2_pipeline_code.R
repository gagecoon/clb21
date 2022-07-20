library(dada2); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

path <- "FASTQ"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq 
# and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), trimRight=c(5,10),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 231:275]

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) ## bad samples = 10-12 CLB2, 12-14 CLB2, 16-18 CLB2, 30-32 CLB2, 32-34 CLB1, 4-6 CLB2, JBpcrBlank

section.size <- 1000     ##this could be any size, so long as it works for your computer. and make sure to alter the refFasta file to be your own file :)

sections <- split(c(1:nrow(seqtab.nochim)),
                  sort(c(1:nrow(seqtab.nochim))%%ceiling(nrow(seqtab.nochim)/section.size)))

sections.taxatest <- lapply(sections,function(x){return(assignTaxonomy(seqtab.nochim[x,],refFasta="~/Desktop/June_2021_CLB_Sediments/16S Processing/silva_nr99_v138.1_train_set.fa",
                                                                       verbose = TRUE))})
taxa.test <- do.call(rbind, sections.taxatest)

is.chloro <- taxa.test[,"Order"] %in% "Chloroplast"
is.mito <- taxa.test[,"Family"] %in% "Mitochondria"

print(dim(seqtab.nochim))
seqtab.nochloro.nomito <- seqtab.nochim[, !(is.chloro | is.mito)]
print(dim(seqtab.nochloro.nomito))

print(dim(taxa.test))
taxa.nochloro.nomito <- taxa.test[!(is.chloro | is.mito), ]
print(dim(taxa.nochloro.nomito))

save(taxa.nochloro.nomito, seqtab.nochloro.nomito,
     file = "CLB_Data_Analysis.rda")

asv.seqs <- colnames(seqtab.nochloro.nomito)
asv_headers <- vector(dim(seqtab.nochloro.nomito)[2], mode = "character")

for (i in 1:dim(seqtab.nochloro.nomito)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep ="_")
}

#files

asv_fasta <- c(rbind(asv_headers, asv.seqs))
write(asv_fasta, "ASV_noMC_CLB21.fa")

asv_tab <- t(seqtab.nochloro.nomito)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASV_counts_noMC_CLB21.csv", sep = ",", quote=F, col.names = NA)

taxa.print <- taxa.nochloro.nomito
rownames(taxa.print) <- NULL
head(taxa.print)

row.names(taxa.print) <- sub(">", "", asv_headers)
write.table(taxa.print, "ASV_tax_noMC_CLB21.csv", sep = ",", quote=F, col.names = NA)



