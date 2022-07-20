library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
library(tidyverse)

theme_set(theme_bw())

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
  labs(x = "Water Content (%)", y = "Depth (cm)") +
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
  labs(x = "Water Content (%)", y = "Depth (cm)") +
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
right = text_grob("Core 1                                                                           Core 2    ", rot = 270))

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
