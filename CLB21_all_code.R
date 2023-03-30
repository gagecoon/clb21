#Gage R. Coon
#gcoon@vols.utk.edu
#gagecoon4c@gmail.com

library(dada2)
library(phyloseq)
library(Biostrings)
library(microbiome)
library(RColorBrewer)
library(rioja)
library(microbiomeutilities)
library(igraph)
library(vegan)
library(stringr)
library(ggpubr)
library(metagMisc)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(knitr)
library(phylosmith)
library(threejs)
library(RCy3)
library(phylosmith)
library(readxl)
library(scales)
library(tidyr)
library(dplyr)
library(xlsx)
library(data.table)
library(circlize)
library(ggalluvial)
library(cluster)
library(factoextra)

theme_set(theme_classic())

#
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

section.size <- 1000     ##this could be any size, so long as it works for your computer

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

### end dada2 pipeline code

##################################
#       Making Data Tables       #
##################################

CLB_ASVs <- read.table(file = 'ASV_counts_noMC_CLB21.csv', sep = ',', 
                       header = TRUE)
colnames(CLB_ASVs) <- sub("X", "", colnames(CLB_ASVs))
row.names(CLB_ASVs) <- CLB_ASVs[,1]

CLB_ASVs <- CLB_ASVs[,-1]  #remove names
CLB_ASVs <- CLB_ASVs[,-43] #remove blank
CLB_ASVs <- CLB_ASVs[,-36] #remove 4-6 CLB2
CLB_ASVs <- CLB_ASVs[,-27] #remove 32-34 CLB1
CLB_ASVs <- CLB_ASVs[,-26] #remove 30-32 CLB2
CLB_ASVs <- CLB_ASVs[,-10] #remove 16-18 CLB2
CLB_ASVs <- CLB_ASVs[,-6]  #remove 12-14 CLB2
CLB_ASVs <- CLB_ASVs[,-4]  #remove 10-12 CLB2

CLB_taxa <- read.table(file = 'ASV_tax_noMC_CLB21.csv', sep = ',', 
                       header = TRUE)
CLB_tax_row <- CLB_taxa[,1]
row.names(CLB_taxa) <- CLB_taxa[,1]
CLB_taxa <- CLB_taxa[,-1]

CLB_sample <- read.table(file = 'CLB21_sample_file.csv', sep = ',', 
                         header = TRUE)
row.names(CLB_sample) <- CLB_sample[,1]
CLB_sample <- CLB_sample[,-1]
CLB_sample$Core <- as.character(CLB_sample$Core)

##################################
#     Making Phyloseq Object     #
##################################

CLB_ASV_table <- otu_table(CLB_ASVs, taxa_are_row = TRUE)
CLB_taxa_table <- tax_table(CLB_taxa)
row.names(CLB_taxa_table) <- CLB_tax_row
colnames(CLB_taxa_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                              "Genus")
CLB_sample_table <- sample_data(CLB_sample)

CLB_ps <- phyloseq(CLB_ASV_table, CLB_sample_table, CLB_taxa_table)

CLB_ps = filter_taxa(CLB_ps, function(x) sum(x) > 0, TRUE)
CLB_ps <- prune_taxa(taxa_sums(CLB_ps) > 5, CLB_ps) #removing those with less than 5 reads
##bac_data code taken from https://github.com/dgiovannelli/SubductCR_16S-diversity/blob/22ae60a5e6d751af7667a04064042cf2f70096fb/Fullerton_et_al_BMS_16S_final_analysis.r
bac_data <- CLB_ps
#List of potential contaminant genera in subsurface 16S rRNA libraries after Sheik et al. 2018 Frontiers in Microbiology
bac_data <- subset_taxa(bac_data, (Genus != "Acinetobacter") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pseudomonas") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Abiotrophia") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Achromobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Acinetobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Actinobacillus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Arcanobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Arcobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Babesia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bacillus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bartonella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bordetella") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Borrelia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Brodetella") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Brucella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Burkholderia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Campylobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Capnocytophaga") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Chlamydia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Clostridium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Comamonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Corynebacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Coxiella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Cronobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Deinococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Dermatophilus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Ehrlichia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Enterococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Erysipelothrix") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Escherichia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Escherichia/Shigella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Flavobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Francisella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Gardnerella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Granulicatella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Haemophilus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Hafnia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Halomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Helicobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Klebsiella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Kocuria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Lawsonia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Legionella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Leptospira") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Listeria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Merkel_cell") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Micrococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Morganella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Mycobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Mycoplasma") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Neisseria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Nocardia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pasteurella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Photobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Plesiomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Propionibacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Proteus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Providencia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pseudomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rhodococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rickettsiae") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Roseomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rothia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Salmonella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Serratia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Shewanella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Shigella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Sphaerophorus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Staphylococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Stenotrophomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Streptococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Treponema") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Vibrio") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Yersinia") | is.na(Genus))

bac_data <- subset_taxa(bac_data,  (Kingdom != "Eukaryota") | is.na(Kingdom))

CLB_ps <- bac_data

CLB_kingdoms <- aggregate_taxa(CLB_ps, "Kingdom")
CLB_kingdoms2 <- subset_taxa(CLB_kingdoms, Kingdom == "Archaea" | Kingdom == "Bacteria")
taxa_sums(CLB_kingdoms2) #archaea vs bacteria percentages

summarize_phyloseq(CLB_ps)

CLB_ps.ra <- transform_sample_counts(CLB_ps, function(x){x / sum(x)}) #relative abundance transformation

###################
# Alpha Diversity #
###################

alpha_meas <- c("Shannon", "Simpson", "Chao1")
richness <- plot_richness(CLB_ps, x="Depth", measures=alpha_meas)
plot(richness)
div <- estimate_richness(CLB_ps, split = FALSE, measures=alpha_meas)
div
CLB_ps1 <- subset_samples(CLB_ps, Core == "1")
CLB_ps2 <- subset_samples(CLB_ps, Core == "2")
div <- estimate_richness(CLB_ps1, split = FALSE, measures=alpha_meas)
div
div <- estimate_richness(CLB_ps2, split = FALSE, measures=alpha_meas)
div

###################
#  Beta Diversity #
###################

ps.prop <- transform_sample_counts(CLB_ps, function(otu) otu/sum(otu))

################
#  Ordination  #
################

ord.pcoa.bray <- ordinate(CLB_ps, method="PCoA", distance = "bray")
ord.nmds.bray <- ordinate(CLB_ps, method="NMDS", distance = "bray")
ord.cca.bray <- ordinate(CLB_ps, method="CCA", distance = "bray")

###################
#       PCoA      #
###################

pcoa <- plot_ordination(CLB_ps, ord.pcoa.bray, color ="Depth", title = "PCoA")
pcoa

###################
#       NMDS      #
###################

##mostly useless code from trying things/analyzing data. 
wh0 = genefilter_sample(CLB_ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(CLB_ps))
CLB_ps.filtered = prune_taxa(wh0, CLB_ps)
CLB_ps.filtered = transform_sample_counts(CLB_ps.filtered, function(x){x / sum(x)})
phylum.sum = tapply(taxa_sums(CLB_ps.filtered), tax_table(CLB_ps.filtered)[, "Phylum"], sum, na.rm=TRUE)
top27phyla = names(sort(phylum.sum, TRUE))[1:27]
CLB_ps.filtered = prune_taxa((tax_table(CLB_ps.filtered)[, "Phylum"] %in% top27phyla), CLB_ps.filtered)
ord.nmds.bray.taxa <- ordinate(CLB_ps.filtered, method="NMDS", distance = "bray")
nmds <- plot_ordination(CLB_ps.filtered, ord.nmds.bray.taxa, type = "taxa", color ="Phylum", title = "NMDS")
nmds
#below is the basic nmds plot in paper
nmds_asv <- plot_ordination(CLB_ps, ord.nmds.bray, color ="Depth", title = "NMDS", shape = "Core") + 
  scale_colour_gradient2(high = "black", mid = "#2458d1", low = "#03cafc", midpoint = 25)
nmds_asv

#checking diversity in the two zones, cutoffs are slightly off
dfmz <- subset_samples(CLB_SRB_MG_RA, Depth >= 34)
dfsmtz <- subset_samples(CLB_SRB_MG_RA, Depth < 34)
sample_data(dfmz)$Zone <- as.character("MZ")
sample_data(dfsmtz)$Zone <- as.character("SMTZ")
dfch4srbzone <- merge_phyloseq(dfmz, dfsmtz)
ord.ch4andsrb <- ordinate(dfch4srbzone, method="NMDS", distance = "bray")
ord.mz <- ordinate(dfmz, method="NMDS", distance = "bray")
ord.smtz <- ordinate(dfsmtz, method="NMDS", distance = "bray")
nmds_mz <- plot_ordination(dfmz, ord.mz, type = "taxa", color ="Genus", title = "NMDS of mz")
nmds_smtz <- plot_ordination(dfsmtz, ord.smtz, type = "taxa", color ="Genus", title = "NMDS of smtz")
nmdszone <- ggarrange(nmds_mz, nmds_smtz)
nmdszone

###################
#       CCA       #
###################

cca <- plot_ordination(CLB_ps, ord.cca.bray, color ="Depth", title = "CCA")
cca

beta_diversity <- ggarrange(pcoa, nmds, cca)
beta_diversity #final plot

##################################
#     New Bar Plot Function      #
##################################

#function to get rid of the plot_bar black lines that show up between asvs. Just an edit in the base code

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

##################################
#  Depth series diversity plots  #
##################################

### Instead of stacked bar plots, here I am making relative abundance plots on the y-axis while
### the x-axis is depth. There should be a different curve with connected dotted lines between 
### each relative abundance point per type of microbe in each group (methanogens vs srb)

### subsetting based on core as the srz, mz, and smtz are not identical between each -
### later these will be separate plots in one figure

CLB_ps1.ra <- subset_samples(CLB_ps.ra, Core == "1") #ps = phyloseq
CLB_ps2.ra <- subset_samples(CLB_ps.ra, Core == "2")

CLB_mg1.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanocellales") 
CLB_mg2.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanomassiliicoccales") 
CLB_mg3.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanomicrobiales") 
CLB_mg4.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanosarciniales")
CLB_mg5.ra1 <- subset_taxa(CLB_ps1.ra, Class == "ANME-1")
CLB_mg6.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanofastidiosales")
CLB_mg_taxa.ra1 <- merge_phyloseq(CLB_mg1.ra1,CLB_mg2.ra1,CLB_mg3.ra1,
                                  CLB_mg4.ra1,CLB_mg5.ra1,CLB_mg6.ra1)

CLB_mg1.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanocellales") 
CLB_mg2.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanomassiliicoccales") 
CLB_mg3.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanomicrobiales") 
CLB_mg4.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanosarciniales")
CLB_mg5.ra2 <- subset_taxa(CLB_ps2.ra, Class == "ANME-1")
CLB_mg6.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanofastidiosales")
CLB_mg_taxa.ra2 <- merge_phyloseq(CLB_mg1.ra2,CLB_mg2.ra2,CLB_mg3.ra2,
                                  CLB_mg4.ra2,CLB_mg5.ra2,CLB_mg6.ra2)
seep1.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "SEEP-SRB1")
seep2.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "SEEP-SRB2")
seep4.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "SEEP-SRB4")
msbl7.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "MSBL7")
sva0081.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "Sva0081 sediment group")
desu.ra1 <- subset_taxa(CLB_ps1.ra, Genus == "Desulfatiglans")
CLB_srb_taxa.ra1 <- merge_phyloseq(seep1.ra1,seep2.ra1,seep4.ra1,msbl7.ra1,sva0081.ra1,desu.ra1)
seep1.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "SEEP-SRB1")
seep2.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "SEEP-SRB2")
seep4.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "SEEP-SRB4")
msbl7.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "MSBL7")
sva0081.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "Sva0081 sediment group")
desu.ra2 <- subset_taxa(CLB_ps2.ra, Genus == "Desulfatiglans")
CLB_srb_taxa.ra2 <- merge_phyloseq(seep1.ra2,seep2.ra2,seep4.ra2,msbl7.ra2,sva0081.ra2,desu.ra2)
### merging otus that are the same order
CLB_mg_taxa.ra1_glom <- tax_glom(CLB_mg_taxa.ra1, "Order")
CLB_mg_taxa.ra2_glom <- tax_glom(CLB_mg_taxa.ra2, "Order")
CLB_srb_taxa.ra1_glom <- tax_glom(CLB_srb_taxa.ra1, "Genus")
CLB_srb_taxa.ra2_glom <- tax_glom(CLB_srb_taxa.ra2, "Genus")
### exporting abundances
write.table(CLB_mg_taxa.ra1_glom %>% psmelt() %>% 
              select(Order, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.order_core1.tsv", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(CLB_mg_taxa.ra2_glom %>% psmelt() %>% 
              select(Order, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.order_core2.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(CLB_srb_taxa.ra1_glom %>% psmelt() %>% 
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.genus_core1.tsv", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(CLB_srb_taxa.ra2_glom %>% psmelt() %>% 
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "ps.relative_abundance.genus_core2.tsv", sep = "\t", quote = F, row.names = F, col.names = T) #exported and 
#manually changed names to be avg depth only so I can plot this now. Made into xlsx file too
ps.relative_abundance.order_core1 <- read_excel("ps.relative_abundance.order_core1.xlsx") #mg 1
ps.relative_abundance.order_core2 <- read_excel("ps.relative_abundance.order_core2.xlsx") #mg 2
ps.relative_abundance.genus_core1 <- read_excel("ps.relative_abundance.genus_core1.xlsx") #srb 1
ps.relative_abundance.genus_core2 <- read_excel("ps.relative_abundance.genus_core2.xlsx") #srb 2
### need to elongate the data so relative abundance is its own column
ps.relative_abundance.order_core1_pivot <-
  pivot_longer(ps.relative_abundance.order_core1, cols = !Order, names_to = "Depth", values_to = "relative_abundance")
ps.relative_abundance.order_core2_pivot <-
  pivot_longer(ps.relative_abundance.order_core2, cols = !Order, names_to = "Depth", values_to = "relative_abundance")
ps.relative_abundance.genus_core1_pivot <-
  pivot_longer(ps.relative_abundance.genus_core1, cols = !Genus, names_to = "Depth", values_to = "relative_abundance")
ps.relative_abundance.genus_core2_pivot <-
  pivot_longer(ps.relative_abundance.genus_core2, cols = !Genus, names_to = "Depth", values_to = "relative_abundance")
### making them as continuous numbers to troubleshoot plot
ps.relative_abundance.order_core1_pivot$Depth <- as.numeric(ps.relative_abundance.order_core1_pivot$Depth)
ps.relative_abundance.order_core1_pivot <- arrange(ps.relative_abundance.order_core1_pivot, desc(Depth))
ps.relative_abundance.order_core1_pivot$relative_abundance <- as.numeric(ps.relative_abundance.order_core1_pivot$relative_abundance)
ps.relative_abundance.order_core2_pivot$Depth <- as.numeric(ps.relative_abundance.order_core2_pivot$Depth)
ps.relative_abundance.order_core2_pivot <- arrange(ps.relative_abundance.order_core2_pivot, desc(Depth))
ps.relative_abundance.order_core2_pivot$relative_abundance <- as.numeric(ps.relative_abundance.order_core2_pivot$relative_abundance)
ps.relative_abundance.genus_core1_pivot$Depth <- as.numeric(ps.relative_abundance.genus_core1_pivot$Depth)
ps.relative_abundance.genus_core1_pivot <- arrange(ps.relative_abundance.genus_core1_pivot, desc(Depth))
ps.relative_abundance.genus_core1_pivot$relative_abundance <- as.numeric(ps.relative_abundance.genus_core1_pivot$relative_abundance)
ps.relative_abundance.genus_core2_pivot$Depth <- as.numeric(ps.relative_abundance.genus_core2_pivot$Depth)
ps.relative_abundance.genus_core2_pivot <- arrange(ps.relative_abundance.genus_core2_pivot, desc(Depth))
ps.relative_abundance.genus_core2_pivot$relative_abundance <- as.numeric(ps.relative_abundance.genus_core2_pivot$relative_abundance)
### plotting
### colors
srb_colorlist = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
### single plots for ggarrange
CLB_mg_lp_core_1 <- ggplot(data = ps.relative_abundance.order_core1_pivot, aes(x = relative_abundance, y = Depth, group = Order)) +
  geom_path(aes(color = Order), linetype = "dashed") +
  geom_point(aes(color = Order, shape = Order)) +
  labs(y = "Depth (cm)", x = "Relative abundance") +
  scale_y_reverse() +
  xlim(0,0.02)
CLB_mg_lp_core_1
CLB_mg_lp_core_2 <- ggplot(data = ps.relative_abundance.order_core2_pivot, aes(x = relative_abundance, y = Depth, group = Order)) +
  geom_path(aes(color = Order), linetype = "dashed") +
  geom_point(aes(color = Order, shape = Order)) +
  labs(y = "Depth (cm)", x = "Relative abundance") +
  scale_y_reverse() +
  xlim(0,0.02)
CLB_mg_lp_core_2
CLB_srb_lp_core_1 <- ggplot(data = ps.relative_abundance.genus_core1_pivot, aes(x = relative_abundance, y = Depth, group = Genus)) +
  geom_path(aes(color = Genus), linetype = "dashed") +
  geom_point(aes(color = Genus, shape = Genus)) +
  labs(y = "Depth (cm)", x = "Relative abundance") +
  xlim(0,0.057) +
  scale_y_reverse() +
  scale_color_manual(values = srb_colorlist)
CLB_srb_lp_core_1
CLB_srb_lp_core_2 <- ggplot(data = ps.relative_abundance.genus_core2_pivot, aes(x = relative_abundance, y = Depth, group = Genus)) +
  geom_path(aes(color = Genus), linetype = "dashed") +
  geom_point(aes(color = Genus, shape = Genus)) +
  labs(y = "Depth (cm)", x = "Relative abundance") +
  scale_y_reverse() +
  xlim(0,0.057) +
  scale_color_manual(values = srb_colorlist)
CLB_srb_lp_core_2
### taking the legends to add to plot
legend_1 <- get_legend(CLB_mg_lp_core_1)
legend_2 <- get_legend(CLB_srb_lp_core_1)
### annotate core number on top graph
CLB_mg_lp_core_1_ann <- annotate_figure(CLB_mg_lp_core_1 + theme(legend.position = "None") + labs(x = NULL), top = text_grob("Core 1"))
CLB_mg_lp_core_2_ann <- annotate_figure(CLB_mg_lp_core_2 + theme(legend.position = "None") + labs(y = NULL, x = NULL), top = text_grob("Core 2"))
### ggarrange
CLB_line_plots <- ggarrange(CLB_mg_lp_core_1_ann, CLB_mg_lp_core_2_ann, 
                            legend_1, CLB_srb_lp_core_1 + theme(legend.position = "None"), 
                            CLB_srb_lp_core_2 + theme(legend.position = "None") + labs(y = NULL), legend_2, ncol = 3, nrow = 2, labels = c("A", "B", "", "C", "D", ""), hjust = 0.01)
CLB_line_plots #final plot, saved as a 5 in by 8 in pdf

##################################
#          Methanogens           #
##################################

CLB_ps1.ra <- subset_samples(CLB_ps.ra, Core == "1")
CLB_ps2.ra <- subset_samples(CLB_ps.ra, Core == "2")

#subsetting various groups into their own variable
anme3 <- subset_taxa(CLB_ps.ra, Genus == "ANME-3")
anme2a <- subset_taxa(CLB_ps.ra, Family=="ANME-2a-2b")
anme2c <- subset_taxa(CLB_ps.ra, Family=="ANME-2c")
anme2a2b2c3 <- merge_phyloseq(anme3,anme2a,anme2c)

#plot of anmes besides anme-1
anme3plot <- plot_bar(anme2a2b2c3, "Depth", fill = "Family", facet_grid = "Core")
anme3plot + coord_flip() + scale_x_reverse()

CLB_mg1.ra <- subset_taxa(CLB_ps.ra, Order == "Methanocellales") 
CLB_mg2.ra <- subset_taxa(CLB_ps.ra, Order == "Methanomassiliicoccales") 
CLB_mg3.ra <- subset_taxa(CLB_ps.ra, Order == "Methanomicrobiales") 
CLB_mg4.ra <- subset_taxa(CLB_ps.ra, Order == "Methanosarciniales")
CLB_mg5.ra <- subset_taxa(CLB_ps.ra, Class == "ANME-1")
CLB_mg6.ra <- subset_taxa(CLB_ps.ra, Order == "Methanofastidiosales")
CLB_mg_taxa.ra <- merge_phyloseq(CLB_mg1.ra,CLB_mg2.ra,CLB_mg3.ra,
                                 CLB_mg4.ra,CLB_mg5.ra,CLB_mg6.ra)
CLB_mg_bar.ra <- my_plot_bar(CLB_mg_taxa.ra,"Depth", fill = "Order", facet_grid = "Core") + coord_flip() + scale_x_reverse() +
  theme(text = element_text(size = 10))
CLB_mg_bar.ra #plot of some methanogens 

CLB_mg1.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanocellales") 
CLB_mg2.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanomassiliicoccales") 
CLB_mg3.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanomicrobiales") 
CLB_mg4.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanosarciniales")
CLB_mg5.ra1 <- subset_taxa(CLB_ps1.ra, Class == "ANME-1")
CLB_mg6.ra1 <- subset_taxa(CLB_ps1.ra, Order == "Methanofastidiosales")
CLB_mg_taxa.ra1 <- merge_phyloseq(CLB_mg1.ra1,CLB_mg2.ra1,CLB_mg3.ra1,
                                  CLB_mg4.ra1,CLB_mg5.ra1,CLB_mg6.ra1)
CLB_mg_bar.ra1 <- my_plot_bar(CLB_mg_taxa.ra1,"Depth", fill = "Order") 
CLB_mg_bar.ra1 + coord_flip() + scale_x_reverse() #just some of core 1 methanogens

CLB_mg1.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanocellales") 
CLB_mg2.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanomassiliicoccales") 
CLB_mg3.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanomicrobiales") 
CLB_mg4.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanosarciniales")
CLB_mg5.ra2 <- subset_taxa(CLB_ps2.ra, Class == "ANME-1")
CLB_mg6.ra2 <- subset_taxa(CLB_ps2.ra, Order == "Methanofastidiosales")
CLB_mg_taxa.ra2 <- merge_phyloseq(CLB_mg1.ra2,CLB_mg2.ra2,CLB_mg3.ra2,
                                  CLB_mg4.ra2,CLB_mg5.ra2,CLB_mg6.ra2)
CLB_mg_bar.ra2 <- my_plot_bar(CLB_mg_taxa.ra2,"Depth", fill = "Order") 
CLB_mg_bar.ra2 + coord_flip() + scale_x_reverse() #just some of core 2 methanogens

m1 <- subset_taxa(CLB_ps, Class == "ANME-1")
m2 <- subset_taxa(CLB_ps, Family == "ANME-2a-2b")
m3 <- subset_taxa(CLB_ps, Family == "ANME-2c")
m4 <- subset_taxa(CLB_ps, Genus == "ANME-3")
m5 <- subset_taxa(CLB_ps, Order == "Methanocellales") 
m6 <- subset_taxa(CLB_ps, Order == "Methanomassiliicoccales") 
m7 <- subset_taxa(CLB_ps, Order == "Methanomicrobiales") 
m8 <- subset_taxa(CLB_ps, Order == "Methanosarciniales")
m9 <- subset_taxa(CLB_ps, Order == "Methanofastidiosales")

ch4anaerobes <- merge_phyloseq(m1,m2,m3,m4,m5,m6,m7,m8,m9) #all methanogens and anme archaea

##################################
#       Complete 16s Data        #
##################################

#giant supplemental plot wiht all the phyla pictured
all_phylum <- tax_glom(CLB_ps, taxrank = "Phylum")
all_phylum.ra <- transform_sample_counts(all_phylum, function(x){x / sum(x)})
data_all_phylum <- psmelt(all_phylum.ra)
data_all_phylum$Phylum <- as.character(data_all_phylum$Phylum)
data_all_phylum$Phylum[data_all_phylum$Abundance < 0.01] <- "< 1% Abundance"
colorCount = length(unique(data_all_phylum$Phylum))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
colorlist_unmodified <- getPalette(colorCount)
colorlist = c("#666666", "#B53446", "#864F70", "#586A9A", "#3881AF", "#3E8E91", "#449C74", "#4AA956", "#58A057", "#6C856F",
              "#806B87", "#95519F", "#AF597D", "#CB6651", "#E77325", "#FF8301", "#FFA60F", "#FFC81D", "#FFEB2B", 
              "#F4EB31", "#DCBD2E", "#C4902B", "#AC6228", "#B55E45", "#CB696D", "#E17596", "#F781BF")
all_tax_plot <- ggplot(data_all_phylum, aes(x=Depth, y=Abundance, fill=Phylum))
all_tax_plot + geom_bar(aes(), stat="identity", position="stack") + facet_grid("Core") +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=15)) + coord_flip() + scale_x_reverse() +
  scale_fill_manual(values = colorlist)

##################################
#        Sulfate-Reducers        #
##################################
#many various groups of srb and ways to visualize. code for graph in paper is further down
seep1 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB1")
seep2 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB2")
seep4 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB4")
msbl7 <- subset_taxa(CLB_ps, Genus == "MSBL7")
sva0081 <- subset_taxa(CLB_ps, Genus == "Sva0081 sediment group")
desu <- subset_taxa(CLB_ps, Genus == "Desulfatiglans")
CLB_SRB <- merge_phyloseq(seep1,seep2,seep4,msbl7,sva0081,desu)
seep1.ra <- subset_taxa(CLB_ps.ra, Genus == "SEEP-SRB1")
seep2.ra <- subset_taxa(CLB_ps.ra, Genus == "SEEP-SRB2")
seep4.ra <- subset_taxa(CLB_ps.ra, Genus == "SEEP-SRB4")
msbl7.ra <- subset_taxa(CLB_ps.ra, Genus == "MSBL7")
sva0081.ra <- subset_taxa(CLB_ps.ra, Genus == "Sva0081 sediment group")
desu.ra <- subset_taxa(CLB_ps.ra, Genus == "Desulfatiglans")
CLB_SRB.ra <- merge_phyloseq(seep1.ra,seep2.ra,seep4.ra,msbl7.ra,sva0081.ra,desu.ra)

CLB_des.ra <-  subset_taxa(CLB_ps.ra, Genus == "Desulfatiglans") 
CLB_see.ra <-  subset_taxa(CLB_ps.ra, Genus == "SEEP-SRB1") 
CLB_sr_taxa.ra <- merge_phyloseq(CLB_des.ra, CLB_see.ra)

CLB_des.ra1 <-  subset_taxa(CLB_ps1.ra, Genus == "Desulfatiglans") 
CLB_see.ra1 <-  subset_taxa(CLB_ps1.ra, Genus == "SEEP-SRB1") 
CLB_sr_taxa.ra1 <- merge_phyloseq(CLB_des.ra1, CLB_see.ra1)

CLB_des.ra2 <-  subset_taxa(CLB_ps2.ra, Genus == "Desulfatiglans") 
CLB_see.ra2 <-  subset_taxa(CLB_ps2.ra, Genus == "SEEP-SRB1") 
CLB_sr_taxa.ra2 <- merge_phyloseq(CLB_des.ra2, CLB_see.ra2)

CLB_sr_bar.ra1 <- plot_bar(CLB_sr_taxa.ra1, "Depth", fill = "Genus",
                           facet_grid = "Genus") 
CLB_sr_bar.ra2 <- plot_bar(CLB_sr_taxa.ra2, "Depth", fill = "Genus",
                           facet_grid = "Genus") 

srgenus <- psmelt(CLB_SRB.ra)
srgenus$Genus <- as.character(srgenus$Genus)
colorCount2 = length(unique(srgenus$Genus))
getPalette2 = colorRampPalette(brewer.pal(6, "Dark2"))
colorlist_unmodified2 <- getPalette2(colorCount2)
colorlist_unmodified2
colorlist2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
CLB_srplot <- my_plot_bar(CLB_SRB.ra, "Depth", fill = "Genus", facet_grid = "Core") + coord_flip() + scale_x_reverse() +
  theme(text = element_text(size = 10)) +  scale_fill_manual(values = colorlist2)
CLB_srplot

srmet <- ggarrange(CLB_mg_bar.ra, CLB_srplot, legend = "bottom", labels = c("A", "B"))
srmet #final paper plot
##################################
#        mg and mt Graphs        #
##################################
#mg = methanogens & mt = methylotrophs

legendplotbasiclegend_mg <- CLB_mg_bar.ra1
legendplot_mg <- legendplotbasiclegend_mg + scale_fill_manual(values= speciesPalette_mg)

legendplotbasiclegend_mt <- CLB_mt_bar.ra1
legendplot_mt <- legendplotbasiclegend_mt + scale_fill_manual(values= speciesPalette_mt)

################################## Methanogens

CLB_mg1.ra <- CLB_mg_bar.ra1 + coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + theme(legend.position = "None") + xlab("Depth (cm)") +
  scale_fill_manual(values= speciesPalette_mg)
CLB_mg2.ra <- CLB_mg_bar.ra2 + coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + rremove("y.text") + xlab("Depth (cm)") +
  theme(legend.position = "None") +
  scale_fill_manual(values= speciesPalette_mg)

################################## Methanotrophs

CLB_mt1.ra <- CLB_mt_bar.ra1 + coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + theme(legend.position = "None") + xlab("Depth (cm)") +
  scale_fill_manual(values= speciesPalette_mt) + ylim(0,0.006)
CLB_mt2.ra <- CLB_mt_bar.ra2 + coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + rremove("y.text") + xlab("Depth (cm)") +
  theme(legend.position = "None") +
  scale_fill_manual(values= speciesPalette_mt) + ylim(0,0.006)

CLB_sr.ra1 <- CLB_sr_bar.ra1 + theme(legend.position = "None") + 
  coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + ylim(0,0.06) + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) + xlab("Depth (cm)")
CLB_sr.ra2 <- CLB_sr_bar.ra2 + theme(legend.position = "None") + 
  coord_flip() + ylab("Relative Abundance") +
  scale_x_reverse() + ylim(0,0.06) + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) + xlab("Depth (cm)")

##################################
#         Universal Plot         #
##################################

g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend_mg <- g_legend(legendplot_mg)
legend_mt <- g_legend(legendplot_mt)

##################################
#            Counts              #
##################################

#counts of groups used in paper

m1 <- subset_taxa(CLB_ps, Class == "ANME-1")
m2 <- subset_taxa(CLB_ps, Family == "ANME-2a-2b")
m3 <- subset_taxa(CLB_ps, Family == "ANME-2c")
m4 <- subset_taxa(CLB_ps, Genus == "ANME-3")
m5 <- subset_taxa(CLB_ps, Order == "Methanocellales") 
m6 <- subset_taxa(CLB_ps, Order == "Methanomassiliicoccales") 
m7 <- subset_taxa(CLB_ps, Order == "Methanomicrobiales") 
m8 <- subset_taxa(CLB_ps, Order == "Methanosarciniales")
m9 <- subset_taxa(CLB_ps, Order == "Methanofastidiosales")

sum(sample_sums(m1)) #12533
sum(sample_sums(m2)) #10
sum(sample_sums(m3)) #67
sum(sample_sums(m4)) #674 total anme = 13289
sum(sample_sums(m5)) #13 = <0.1%
sum(sample_sums(m6)) #616 = 2.6%
sum(sample_sums(m7)) #5993 = 25.4%
sum(sample_sums(m8)) #7515 - ANME @ 751 = 6764 = 28.7%
sum(sample_sums(m9)) #10189 total methanogens - anme = 23575 = 43.2%

seep1 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB1")
seep2 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB2")
seep4 <- subset_taxa(CLB_ps, Genus == "SEEP-SRB4")
msbl7 <- subset_taxa(CLB_ps, Genus == "MSBL7")
sva0081 <- subset_taxa(CLB_ps, Genus == "Sva0081 sediment group")
desu <- subset_taxa(CLB_ps, Genus == "Desulfatiglans")

sum(sample_sums(seep1))
sum(sample_sums(seep2))
sum(sample_sums(seep4))
sum(sample_sums(msbl7))
sum(sample_sums(sva0081))
sum(sample_sums(desu))

srbmet <- merge_phyloseq(m1,m2,m3,m4,m5,m6,m7,m8,m9,seep1,seep2,seep4,msbl7,sva0081,desu)
srbmet.ra <- transform_sample_counts(srbmet, function(x){x / sum(x)})
abundance_heatmap(srbmet.ra,classification = "Genus", treatment = "Depth")
variable_correlation_heatmap(srbmet.ra,classification = "Genus",variables = "Depth",method = "spearman")

# network analysis for srb and methanogens with no zone seperation - via cytoscape (not used)

plot_net(srbmet, type = "taxa", maxdist = 0.7, color = "Genus", point_size = 2, point_alpha = .75)
srbmet_network <- make_network(srbmet, type = "taxa", max.dist = 0.6)
plot_network(srbmet_network, srbmet, type = "taxa", color = "Genus")
write_graph(srbmet_network, "/Users/gagercoon/Desktop/network.txt", "gml")
createNetworkFromIgraph(srbmet_network,
                        title = "From igraph",
                        collection = "My Igraph Network Collection")

#heatmap data smtz == 34-42 core 1 and 40-42 core 2 

srbmet1 <- subset_samples(srbmet, Core = 1)
srbmet2 <- subset_samples(srbmet, Core = 2)

#smtz data for network analysis and heatmap

srbmet_smtz1 <- subset_samples(srbmet1, Depth >= 34)
srbmet_smtz2 <- subset_samples(srbmet2, Depth >= 40)
srbmet_smtz <- merge_phyloseq(srbmet_smtz1, srbmet_smtz2)

srbmet_smtz <- filter_taxa(srbmet_smtz, function(x) sum(x) > 0, TRUE)
srbmet_smtz <- prune_taxa(taxa_sums(srbmet_smtz) > 5, srbmet_smtz)
srbmet_smtz_network <- make_network(srbmet_smtz, type = "taxa", max.dist = 0.5)
plot_network(srbmet_smtz_network, srbmet_smtz, type = "taxa", color = "Genus")

write_graph(srbmet_smtz_network, "/Users/gagercoon/Desktop/smtz_network_jul15.txt", "gml")
createNetworkFromIgraph(srbmet_smtz_network,
                        title = "smtz",
                        collection = "CLB21 smtz vs srz")

#heatmap export
smtz_network_otu <- as(otu_table(srbmet_smtz), "matrix")
smtz_network_otu <- as.data.frame(smtz_network_otu)
write.csv(smtz_network_otu, "/Users/gagercoon/Desktop/smtz_network_otu_jul15.csv")
smtz_network_taxa <- as(tax_table(srbmet_smtz), "matrix")
smtz_network_taxa <- as.data.frame(smtz_network_taxa)
write.csv(smtz_network_taxa, "/Users/gagercoon/Desktop/smtz_network_taxa_jul15.csv")
#went in and edited excel file to have names in second row by asv label and counts


#srz data for network analysis and heatmap
srbmet_srz1 <- subset_samples(srbmet1, Depth <= 34)
srbmet_srz2 <- subset_samples(srbmet2, Depth <= 40)
srbmet_srz <- merge_phyloseq(srbmet_srz1, srbmet_srz2)

srbmet_srz <- filter_taxa(srbmet_srz, function(x) sum(x) > 0, TRUE)
srbmet_srz <- prune_taxa(taxa_sums(srbmet_srz) > 5, srbmet_srz)
srbmet_srz_network <- make_network(srbmet_srz, type = "taxa", max.dist = 0.5)
plot_network(srbmet_srz_network, srbmet_srz, type = "taxa", color = "Genus")

write_graph(srbmet_srz_network, "/Users/gagercoon/Desktop/srz_network_jul15.txt", "gml")
createNetworkFromIgraph(srbmet_srz_network,
                        title = "srz",
                        collection = "CLB21 smtz vs srz")
#heatmap export
srz_network_otu <- as(otu_table(srbmet_srz), "matrix")
srz_network_otu <- as.data.frame(srz_network_otu)
write.csv(srz_network_otu, "/Users/gagercoon/Desktop/srz_network_otu_jul15.csv")
srz_network_taxa <- as(tax_table(srbmet_srz), "matrix")
srz_network_taxa <- as.data.frame(srz_network_taxa)
write.csv(srz_network_taxa, "/Users/gagercoon/Desktop/srz_network_taxa_jul15.csv")
#went in and edited excel file to have names in second row by asv label and counts

### EA plots
EA_data <- read_excel("/Users/gagercoon/Desktop/CLB Paper Data/CLB21 EA data.xls", sheet = 2)
BH_char <- as.character(EA_data$BH)
depth <- EA_data$top_depth
total_carbon <- EA_data$"Total Carbon % (untreated)"
inorganic_carbon <- EA_data$"Inorganic Carbon % (untreated - treated)"
organic_carbon <- EA_data$"Organic Carbon (treated)"
total_nitrogen <- EA_data$"Total Nitrogen %"
inorganic_nitrogen <- EA_data$"Nontreated Nitrogen %"
organic_nitrogen <- EA_data$"Treated Nitrogen %"
total_d13C <- EA_data$"d13C Total (untreated)"
inorganic_d13C <- EA_data$"d13C Inorganic (untreated - treated)"
organic_d13C <- EA_data$"d13C Organic (treated)"
total_d15N <- EA_data$"d15N Total (untreated)"
dna_concentrations <- EA_data$"DNA (ng/uL)"
organic_c_n_ratio <- EA_data$"C/N Treated"
inorganic_c_n_ratio <- EA_data$"C/N Nontreated"
porosity <- EA_data$"Water content %"

total_carbon_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = total_carbon), aes(x = total_carbon, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
inorganic_carbon_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = inorganic_carbon), aes(x = inorganic_carbon, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
organic_carbon_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = organic_carbon), aes(x = organic_carbon, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
total_nitrogen_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = total_nitrogen), aes(x = total_nitrogen, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
inorganic_nitrogen_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = inorganic_nitrogen), aes(x = inorganic_nitrogen, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
organic_nitrogen_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = organic_nitrogen), aes(x = organic_nitrogen, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
total_d13C_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = total_d13C), aes(x = total_d13C, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
inorganic_d13C_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = inorganic_d13C), aes(x = inorganic_d13C, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
organic_d13C_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = organic_d13C), aes(x = organic_d13C, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
total_d15N_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = total_d15N), aes(x = total_d15N, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
dna_concentrations_plot <- ggplot(data = EA_data, aes(x = dna_concentrations, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
organic_c_n_ratio_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = organic_c_n_ratio), aes(x = organic_c_n_ratio, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
inorganic_c_n_ratio_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = inorganic_c_n_ratio), aes(x = inorganic_c_n_ratio, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()
porosity_plot <- ggplot(data = remove_missing(EA_data, na.rm = TRUE, vars = porosity), aes(x = porosity, y = depth)) +
  geom_point(aes(color = BH_char)) + scale_y_reverse()


ggarrange(total_carbon_plot, organic_carbon_plot, 
          total_nitrogen_plot, organic_nitrogen_plot,
          total_d13C_plot, organic_d13C_plot, nrow = 3, ncol = 2, common.legend = TRUE)

#start geochem 
#separating all data

geo_data <- read_excel("CLB21_Geochemistry.xlsx", sheet = 1)
depth <- geo_data$Depth
sulfide <- geo_data$`Sulfide (mM)`
porosity <- geo_data$Porosity
porositysd <- geo_data$`Porosity SD`
methane <- geo_data$`Methane (mM)`
methanesd <- geo_data$`Methane SD`
core <- as.character(geo_data$Core)
label_percent()(porosity)
sulfate <- geo_data$`Sulfate (mM)`
dna <- geo_data$DNA

core.data.1 <- read_excel("CLB21_Geochemistry.xlsx", sheet = 2)
depth.1 <- core.data.1$Depth
sulfide.1 <- core.data.1$`Sulfide (mM)`
porosity.1 <- core.data.1$Porosity
porositysd.1 <- core.data.1$`Porosity SD`
methane.1 <- core.data.1$`Methane (mM)`
methanesd.1 <- core.data.1$`Methane SD`
core.1 <- as.character(core.data.1$Core)
porosity.1.p <- label_percent()(porosity.1)
porosity.1.p <- parse_number(porosity.1.p, character())
porositysd.1.p <- label_percent()(porositysd.1)
porositysd.1.p <- parse_number(porositysd.1.p, character())
sulfate.1 <- core.data.1$`Sulfate (mM)`
dna.1 <- core.data.1$DNA

core.data.2 <- read_excel("CLB21_Geochemistry.xlsx", sheet = 3)
depth.2 <- core.data.2$Depth
sulfide.2 <- core.data.2$`Sulfide (mM)`
porosity.2 <- core.data.2$Porosity
porositysd.2 <- core.data.2$`Porosity SD`
methane.2 <- core.data.2$`Methane (mM)`
methanesd.2 <- core.data.2$`Methane SD`
core.2 <- as.character(core.data.2$Core)
porosity.2.p <- label_percent()(porosity.2)
porosity.2.p <- parse_number(porosity.2.p, character())
porositysd.2.p <- label_percent()(porositysd.2)
porositysd.2.p <- parse_number(porositysd.2.p, character())
sulfate.2 <- core.data.2$`Sulfate (mM)`
dna.2 <- core.data.2$DNA

#graphs with core 1 and 2

ggporosity <- ggplot(data = geo_data, aes(x = porosity, y = depth)) +
  geom_point(aes(color = core)) +
  ggtitle("Porosity") +
  labs(x = "Water Content", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  geom_errorbar(aes(xmin = porosity - porositysd, xmax = porosity + porositysd))

ggsulfide <- ggplot(data = geo_data, aes(x = sulfide, y = depth)) +
  geom_point(aes(color = core)) +
  ggtitle("Sulfide") +
  labs(x = "Sulfide (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core")

ggmethane <- ggplot(data = geo_data, aes(x = methane, y = depth)) +
  geom_point(aes(color = core)) +
  ggtitle("Methane") +
  labs(x = "Methane (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  geom_errorbar(aes(xmin = methane - methanesd, xmax = methane + methanesd))

ggsulfate <- ggplot(data = geo_data, aes(x = sulfate, y = depth)) +
  geom_point(aes(color = core)) +
  ggtitle("Sulfate") +
  labs(x = "Sulfate (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core")

ggdna <- ggplot(data = geo_data, aes(x = dna, y = depth)) +
  geom_point(aes(shape = cut(dna, c(-Inf, 0.1, 59.9, Inf)))) +
  ggtitle(" ") +
  labs(x = "DNA (ng/µL)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0,60) +
  scale_shape_manual(values = c("(-Inf,0.1]" = 4,
                                "(0.1,59.9]" = 19,
                                "(59.9, Inf]" = 8),
                     labels = c("too low", "0.1-59.9", "too high"))

#grid.arrange(ggporosity, ggsulfide, ggmethane, ggsulfate, nrow = 2)

#core 1 graphs

ggporosity.1 <- ggplot(data = core.data.1, aes(x = porosity.1.p, y = depth.1)) +
  geom_errorbar(aes(xmin = porosity.1.p - porositysd.1.p, xmax = porosity.1.p + porositysd.1.p),
                color = "#87C9ED") +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Water Content %", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(45, 75)

ggsulfide.1 <- ggplot(data = core.data.1, aes(x = sulfide.1, y = depth.1)) +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Sulfide (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  xlim(0,10.5) +
  labs(color="Core")

ggmethane.1 <- ggplot(data = core.data.1, aes(x = methane.1, y = depth.1)) +
  geom_errorbar(aes(xmin = methane.1 - methanesd.1, xmax = methane.1 + methanesd.1),
                color = "#87C9ED") +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Methane (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0, 1.5)

ggsulfate.1 <- ggplot(data = core.data.1, aes(x = sulfate.1, y = depth.1)) +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Sulfate (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0,22)

ggdna.1 <- ggplot(data = core.data.1, aes(x = dna.1, y = depth.1)) +
  geom_point(aes(shape = cut(dna.1, c(-Inf, 0.1, 59.9, Inf)))) +
  ggtitle(" ") +
  labs(x = "DNA (ng/µL)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0,60) +
  scale_shape_manual(values = c("(-Inf,0.1]" = 4,
                                "(0.1,59.9]" = 19,
                                "(59.9, Inf]" = 8))

#grid.arrange(ggporosity.1, ggsulfide.1, ggmethane.1, ggsulfate.1, nrow = 2)

#core 2 graphs

ggporosity.2 <- ggplot(data = core.data.2, aes(x = porosity.2.p, y = depth.2)) +
  geom_errorbar(aes(xmin = porosity.2.p - porositysd.2.p, xmax = porosity.2.p + porositysd.2.p),
                color = "#87C9ED") +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Water Content %", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(45, 75)

ggsulfide.2 <- ggplot(data = core.data.2, aes(x = sulfide.2, y = depth.2)) +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Sulfide (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  xlim(0,10.5) +
  labs(color="Core")

ggmethane.2 <- ggplot(data = core.data.2, aes(x = methane.2, y = depth.2)) +
  geom_errorbar(aes(xmin = methane.2 - methanesd.2, xmax = methane.2 + methanesd.2),
                color = "#87C9ED") +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Methane (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0, 1.5)

ggsulfate.2 <- ggplot(data = core.data.2, aes(x = sulfate.2, y = depth.2)) +
  geom_point(aes()) +
  ggtitle(" ") +
  labs(x = "Sulfate (mM)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0,22)

ggdna.2 <- ggplot(data = core.data.2, aes(x = dna.2, y = depth.2)) +
  geom_point(aes(shape = cut(dna.2, c(-Inf, 0.1, 59.9, Inf)))) +
  ggtitle(" ") +
  labs(x = "DNA (ng/µL)", y = "Depth (cm)") +
  scale_y_continuous(trans = "reverse") +
  labs(color="Core") +
  xlim(0,60) +
  scale_shape_manual(values = c("(-Inf,0.1]" = 4,
                                "(0.1,59.9]" = 19,
                                "(59.9, Inf]" = 8))

#grid.arrange(ggporosity.2, ggsulfide.2, ggmethane.2, ggsulfate.2, nrow = 2)

#all core 1 and core 2 graphs

geo <- ggarrange(ggporosity.1 + theme(plot.title = element_text(hjust = 0.5),
                                      axis.title.x = element_text(vjust=-0.2)),
                 ggsulfide.1 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_text(vjust=-0.2)),
                 ggmethane.1 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_text(vjust=-0.2)), 
                 ggsulfate.1 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5)),
                 ggdna.1 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none"),
                 ggporosity.2 + theme(plot.title = element_text(hjust = 0.5),
                                      axis.title.x = element_text(vjust=-0.2)),
                 ggsulfide.2 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_text(vjust=-0.2)),
                 ggmethane.2 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_text(vjust=-0.2)),
                 ggsulfate.2 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5)),
                 ggdna.2 + labs(y = NULL) + 
                   theme(axis.text.y=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none"),
                 labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), nrow = 2, ncol = 5)

annotate_figure(geo,
                right = text_grob("Core 1                                                       Core 2    ", rot = 270))

#grid.arrange(ggsulfide.1, ggsulfide.2)

#sulfate and sulfide overlay

ScaleFactor <- max(methane, na.rm = TRUE)/max(sulfate, na.rm = TRUE)

ggsulfurplot <- ggplot(geo_data, aes(y = depth)) + 
  geom_point(aes(x = methane), color = "blue") +
  geom_point(aes(x = sulfate * ScaleFactor), color = "black") + 
  scale_x_continuous(name="Methane (mM)", sec.axis=sec_axis(~./ScaleFactor, name="Sulfate (mM)")) +
  theme(
    axis.title.x.bottom=element_text(color="blue"),
    axis.text.x.bottom=element_text(color="blue"),
    axis.title.x.top=element_text(color="black"),
    axis.text.x.top=element_text(color="black")
  ) +
  scale_y_continuous(trans = "reverse") +
  facet_grid(col=vars(Core))

ggsulfurplot
#correlation
cor.test(sulfate, sulfide, method="pearson", use = "complete.obs")
ggscatter(geo_data, x = "Sulfate (mM)", y = "Sulfide (mM)", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Sulfate", ylab = "Sulfide")
shapiro.test(sulfate)
shapiro.test(sulfide)
ggqqplot(sulfate, ylab = "Sulfate")
ggqqplot(sulfide, ylab = "Sulfide")

#heatmap code
#main code source from Timothy J. Rogers
#df from file output from heatmap export section of main data analysis
CLB_df <- read_csv("/Users/gagercoon/Desktop/Heatmap CSVs July 15th CLB/srz_network_otu_jul15.csv")


CLB_gather_2 <- gather(CLB_df, Site, Abundance, 3:36)
#view(CLB_gather_2)

CLB_select_test <- CLB_gather_2[c("Site","Order","Abundance")]
#view(CLB_select_test)

CLB_sum <- CLB_select_test%>%
  group_by(Site, Order) %>%
  summarise(total = sum(Abundance))
#view(CLB_sum)

CLB_spread <- spread(CLB_sum, Site , total, fill = NA)
CLB_spread[is.na(CLB_spread)] <- 0
#view(CLB_spread)


test_clb <- as.dendrogram(
  hclust(d= dist(cor(CLB_spread[,-1], method = "spearman"))))
#plot(test_clb)

CLB_order <- order.dendrogram(test_clb)
#view(CLB_order)
cor_CLB <- CLB_spread[,c(1,(CLB_order+length(1)))]
#view(cor_CLB)

cor_CLB_2 <- cor_CLB %>% column_to_rownames(., "Order")
cor_CLB_3  <- as.data.frame(t(cor_CLB_2))
#view(cor_CLB_3)

taxon_hc <- as.dendrogram(
  hclust(d= dist(cor(cor_CLB_3, method = "spearman"))))
plot(taxon_hc)

#view(cor(cor_CLB_3, method = "spearman"))
taxon_correlation = as.data.frame(cor(cor_CLB_3, method = "spearman"))
#view(taxon_correlation)

row_name_test <- tibble::rownames_to_column(taxon_correlation, "taxon1")
#view(row_name_test)

taxon_gather <- gather(row_name_test, taxon2, correlation, 2:16)
#view(taxon_gather)

taxon_correlation_test <- ggplot(data=taxon_gather, aes(x=taxon1, y=taxon2, fill=correlation)) +
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(size = 12)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red"), limits = c(-0.8,1)) +
  theme(legend.position = "right") +
  scale_y_discrete(position = "right")

print(taxon_correlation_test)

#smtz

QCLB_df <- read_csv("/Users/gagercoon/Desktop/Heatmap CSVs July 15th CLB/smtz_network_otu_jul15.csv")

QCLB_gather_2 <- gather(QCLB_df, Site, Abundance, 3:10)
#view(CLB_gather_2)

QCLB_select_test <- QCLB_gather_2[c("Site","Order","Abundance")]
#view(CLB_select_test)

QCLB_sum <- QCLB_select_test%>%
  group_by(Site, Order) %>%
  summarise(total = sum(Abundance))
#view(CLB_sum)

QCLB_spread <- spread(QCLB_sum, Site , total, fill = NA)
QCLB_spread[is.na(QCLB_spread)] <- 0
#view(CLB_spread)


Qtest_clb <- as.dendrogram(
  hclust(d= dist(cor(QCLB_spread[,-1], method = "spearman"))))
#plot(Qtest_clb)

QCLB_order <- order.dendrogram(Qtest_clb)
#view(CLB_order)
Qcor_CLB <- QCLB_spread[,c(1,(QCLB_order+length(1)))]
#view(cor_CLB)

Qcor_CLB_2 <- Qcor_CLB %>% column_to_rownames(., "Order")
Qcor_CLB_3  <- as.data.frame(t(Qcor_CLB_2))
#view(cor_CLB_3)

Qtaxon_hc <- as.dendrogram(
  hclust(d= dist(cor(Qcor_CLB_3, method = "spearman"))))
plot(Qtaxon_hc)

#view(cor(cor_CLB_3, method = "spearman"))
Qtaxon_correlation = as.data.frame(cor(Qcor_CLB_3, method = "spearman"))
#view(taxon_correlation)

Qrow_name_test <- tibble::rownames_to_column(Qtaxon_correlation, "taxon1")
#view(row_name_test)

Qtaxon_gather <- gather(Qrow_name_test, taxon2, correlation, 2:13)
#view(taxon_gather)

Qtaxon_correlation_test <- ggplot(data=Qtaxon_gather, aes(x=taxon1, y=taxon2, fill=correlation)) +
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "bottom") +
  theme(text = element_text(size = 12)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red"), limits = c(-0.8,1)) +
  theme(legend.position = "right") +
  scale_y_discrete(position = "right")

print(Qtaxon_correlation_test)

ggarrange(taxon_correlation_test, Qtaxon_correlation_test, common.legend = TRUE) # main plot

