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

theme_set(theme_classic())

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
nmds_asv <- plot_ordination(CLB_ps, ord.nmds.bray, color ="Depth", title = "NMDS")
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
CLB_mg_bar.ra <- my_plot_bar(CLB_mg_taxa.ra,"Depth", fill = "Order", facet_grid = "Core") + coord_flip() + scale_x_reverse()
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

CLB_srplot <- my_plot_bar(CLB_SRB.ra, "Depth", fill = "Genus", facet_grid = "Core") + coord_flip() + scale_x_reverse()
CLB_srplot

srmet <- ggarrange(CLB_mg_bar.ra, CLB_srplot, legend = "bottom")
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

