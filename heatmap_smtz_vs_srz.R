#main code source from Timothy J. Rogers

library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
library(readxl)
library(dplyr)
library(xlsx)
library(data.table)
library(circlize)
library(ggalluvial)
library(cluster)
library(factoextra)

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
