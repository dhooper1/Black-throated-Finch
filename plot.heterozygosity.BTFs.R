## Load packages
library(tidyverse)
library(ggridges)
library(MetBrewer)
library(cowplot)

setwd("/Users/danielhooper/Documents/Projects/BlackthroatedFinch/Manuscript/Analyses/heterozygosity")
rm(list = ls())
theme_set(theme_bw())

###
### Read in heterozygosity results for ddRAD and WGS datasets
###
ddRAD <- read_delim("ddRAD/heterozygosity.ddRAD.miss20.csv", delim = ",", col_names = c("IID", "FID", "TAXON", "POP", "LONGITUDE", "LATITUDE", "SITES", "HET.SITES", "HETEROZYGOSITY"), skip = 1)
ddRAD.pop <- read_delim("ddRAD/heterozygosity.ddRAD.miss20.pop.csv", delim = ",", col_names = c("FID", "TAXON", "POP", "LONGITUDE", "LATITUDE", "N", "MEAN.HET"), skip = 1)
WGS <- read_delim("WGS/heterozygosity.WGS.csv", delim = ",", col_names = c("IID", "FID", "TAXON", "POP", "LONGITUDE", "LATITUDE", "SEX", "HETEROZYGOSITY"), skip = 1)
WGS.pop <- read_delim("WGS/heterozygosity.WGS.pop.csv", delim = ",", col_names = c("FID", "TAXON", "POP", "LONGITUDE", "LATITUDE", "N", "MEAN.HET"), skip = 1)

## Evaluate summaries by taxon
ddRAD |> group_by(TAXON) |> summarize(mean = mean(HETEROZYGOSITY), stdev = sd(HETEROZYGOSITY))
WGS |> group_by(TAXON) |> summarize(mean = mean(HETEROZYGOSITY), stdev = sd(HETEROZYGOSITY))

## MODEL TEST
## Is heterozygosity significantly different between atropygialis, cincta, and hecki?
## ddRAD dataset
model1 <- lm(ddRAD$HETEROZYGOSITY ~ ddRAD$TAXON)
ANOVA1 <- aov(model1)
## Tukey test to study each pair of treatment :
TUKEY1 <- TukeyHSD(x=ANOVA1, 'ddRAD$TAXON', conf.level=0.95)
plot(TUKEY1 , las=1 , col="red")
## WGS dataset
model2 <- lm(WGS$HETEROZYGOSITY ~ WGS$TAXON)
ANOVA2 <- aov(model2)
## Tukey test to study each pair of treatment :
TUKEY2 <- TukeyHSD(x=ANOVA2, 'WGS$TAXON', conf.level=0.95)
plot(TUKEY2 , las=1 , col="red")

## Find common samples between the WGS and ddRAD datasets to test correlation between estimates
common_samples <- intersect(WGS$IID, ddRAD$IID)

WGS.subset <- WGS |> filter(IID %in% common_samples)
WGS.subset <- WGS.subset |> arrange(IID)
ddRAD.subset <- ddRAD |> filter(IID %in% common_samples)
ddRAD.subset <- ddRAD.subset |> arrange(IID)

df <- WGS.subset |> select(IID, TAXON)
df$WGS.het <- WGS.subset$HETEROZYGOSITY
df$ddRAD.het <- ddRAD.subset$HETEROZYGOSITY
df <- df |> filter(IID != "BT51")

df |> ggplot(aes(x = WGS.het, y = ddRAD.het, color = TAXON)) +
  geom_point()

## How well correlated are the estimates of heterozygosity?
cor.test(df$WGS.het, df$ddRAD.het, method=c("pearson")) #-0.1967188
het.model <- lm(df$WGS.het ~ df$ddRAD.het)
summary(het.model)


## Plot population variation as boxplots

## All three taxa: hecki, atropygialis, cincta
taxon_colors <- c("#FAB255", "#DD5129", "#2c7fb8")
ddRAD |> filter(IID != "BT51") |> mutate(TAXON = fct_relevel(TAXON, "hecki", "atropygialis", "cincta")) |> ggplot(aes(x = as.factor(TAXON), y = HETEROZYGOSITY, fill = TAXON, group = TAXON)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) + 
  scale_color_manual(values=taxon_colors) +
  geom_boxplot(aes(alpha = 0.9), outlier.shape = NA) +
  scale_fill_manual(values=taxon_colors) +
  labs(x = "Species", y = "Mean Site-Based Heterozygosity", title = "ddRAD Dataset: 53,049 SNPs") +
  ylim(0.03, 0.07) +
  theme(legend.position = "none")

WGS |> filter(IID != "BT51" & IID != "BT01" & IID != "BT29" & IID != "LH14") |> mutate(TAXON = fct_relevel(TAXON, "hecki", "atropygialis", "cincta")) |> ggplot(aes(x = as.factor(TAXON), y = HETEROZYGOSITY, fill = TAXON, group = TAXON)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) + 
  scale_color_manual(values=taxon_colors) +
  geom_boxplot(aes(alpha = 0.9), outlier.shape = NA) +
  scale_fill_manual(values=taxon_colors) +
  labs(x = "Species", y = "Mean Site-Based Heterozygosity", title = "WGS Dataset") +
  ylim(0.00325, 0.005) +
  theme(legend.position = "none")

## Just cincta populations
pop_colors <- c("#FAB255", "#ABB269", "#5DB27D", "#34A28C", "#9e9ac8", "#2C7590", "#84635C", "#DD5129")
ddRAD |> filter(TAXON == "cincta") |> mutate(FID = fct_relevel(FID, "SC", "AW", "FD", "MD", "DD", "CD", "YD", "AD")) |> ggplot(aes(x = as.factor(FID), y = HETEROZYGOSITY, fill = FID, group = FID)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) + 
  scale_color_manual(values=pop_colors) +
  geom_boxplot(aes(alpha = 0.9), outlier.shape = NA) +
  scale_fill_manual(values=pop_colors) +
  labs(x = "Population", y = "Mean Site-Based Heterozygosity", title = "ddRAD Dataset: 53,049 SNPs") +
  ylim(0.03, 0.07) +
  theme(legend.position = "none")

## Just atropygialis populations [N >= 3]
ddRAD |> filter(TAXON == "atropygialis" & FID != "MR" & FID != "KL" & FID != "HB" & FID != "CG" & FID != "CR") |> mutate(FID = fct_relevel(FID, "HR", "NP", "RL", "PP", "KO", "IS", "ME", "DR")) |> ggplot(aes(x = as.factor(FID), y = HETEROZYGOSITY, fill = FID, group = FID)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) + 
  scale_color_manual(values=pop_colors) +
  geom_boxplot(aes(alpha = 0.9), outlier.shape = NA) +
  scale_fill_manual(values=pop_colors) +
  labs(x = "Population", y = "Mean Site-Based Heterozygosity", title = "ddRAD Dataset: 53,049 SNPs") +
  ylim(0.03, 0.07) +
  theme(legend.position = "none")

## MODEL TEST
## Is heterozygosity significantly different between cincta populations?
## ddRAD dataset
cincta <- ddRAD |> filter(TAXON == "cincta")
model3 <- lm(cincta$HETEROZYGOSITY ~ cincta$FID)
ANOVA3 <- aov(model3)
## Tukey test to study each pair of treatment :
TUKEY3 <- TukeyHSD(x=ANOVA3, 'cincta$FID', conf.level=0.95)
plot(TUKEY3 , las=1 , col="red")

## MODEL TEST
## Is heterozygosity significantly different between atropygialis, cincta, and hecki if you remove western Townsville pops?
## ddRAD dataset
tmp.ddRAD <- ddRAD |> filter(FID != "CD" & FID != "YD")
model <- lm(tmp.ddRAD$HETEROZYGOSITY ~ tmp.ddRAD$TAXON)
ANOVA <- aov(model)
## Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'tmp.ddRAD$TAXON', conf.level=0.95)
plot(TUKEY , las=1 , col="red")
