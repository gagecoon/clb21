library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
library(tidyverse)

theme_set(theme_bw())

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
