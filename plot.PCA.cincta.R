## Load packages
library(tidyverse)
library(broom)
library(MetBrewer)
library(cowplot)

setwd("/Users/danielhooper/Documents/Projects/BlackthroatedFinch/Manuscript/Analyses/PCA")
rm(list = ls())
theme_set(theme_bw())

## Plot PCA results for sequentially more closely related groups of cincta

## Read in PCA results
pca1 <- read.csv("ddRAD.cincta.all.mm85.mac2.1kb_thin.pca1-10.eigenvec.csv", header = TRUE)
pca2 <- read.csv("ddRAD.cincta.minus_GB.mm85.mac2.1kb_thin.pca1-10.eigenvec.csv", header = TRUE)
pca3 <- read.csv("ddRAD.cincta.minus_GB_SC.mm85.mac2.1kb_thin.pca1-10.eigenvec.csv", header = TRUE)

#Select color palettes for optimized figure look
pop_colors = met.brewer(name="Egypt", n=8, type="continuous", direction=-1)
pop_colors1 <- c("#FAB255", "#ABB269", "#5DB27D", "#34A28C", "#9e9ac8", "#2C7590", "#84635C", "#DD5129")
pop_colors2 <- c("#FAB255", "#ABB269", "#5DB27D", "#34A28C", "#9e9ac8", "#2C7590", "#84635C")
pop_colors3 <- c("#ABB269", "#5DB27D", "#34A28C", "#9e9ac8", "#2C7590", "#84635C")

## Panel #1
a <- pca1 |> mutate(FID = fct_relevel(FID, "SC", "AW", "FD", "MD", "DD", "CD", "YD", "GB")) |> ggplot(aes(x=-PC1, y=PC2, colour = FID)) + 
  geom_point(aes(alpha = 0.95), size=5) + 
  scale_color_manual(values=pop_colors1) + 
  labs(x = "PC1 (14.5%)", y = "PC2 (12.0%)") +
  geom_point(shape = 1, size = 5.0, colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ggtitle("All cincta samples", subtitle = "106177 SNPs")

## Panel #2
b <- pca2 |> mutate(FID = fct_relevel(FID, "SC", "AW", "FD", "MD", "DD", "CD", "YD")) |> ggplot(aes(x=PC1, y=-PC2, colour = FID)) + 
  geom_point(aes(alpha = 0.95), size=5) + 
  scale_color_manual(values=pop_colors2) + 
  labs(x = "PC1 (14.0%)", y = "PC2 (11.2%)") +
  geom_point(shape = 1, size =5.0, colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ggtitle("All cincta samples - GB excluded", subtitle = "91319 SNPs")

## Panel #3
c <- pca3 |> mutate(FID = fct_relevel(FID, "AW", "FD", "MD", "DD", "CD", "YD")) |> ggplot(aes(x=PC1, y=PC2, colour = FID)) + 
  geom_point(aes(alpha = 0.95), size=5) + 
  scale_color_manual(values=pop_colors3) + 
  labs(x = "PC1 (14.0%)", y = "PC2 (11.5%)") +
  geom_point(shape = 1, size = 5.0, colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ggtitle("All cincta samples - GB and SC excluded", subtitle = "90509 SNPs")

## Plot all three panels as a single row
plot_grid(a, b, c, nrow = 1)
