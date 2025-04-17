## Load packages
library(tidyverse)
library(ggridges)
library(MetBrewer)
library(cowplot)
library(qqman)

setwd("/Users/danielhooper/Documents/Projects/BlackthroatedFinch/Manuscript/Analyses/Fst")
rm(list = ls())
theme_set(theme_bw())

## Read in vcftools  20kb sliding window - with 10kb step size - output for Fst
order <- c("chr2", "chr1", "chr3", "chr4", "chr1A", "chr5", "chr6","chr7","chr8","chr9","chr10", "chr11", "chr12", "chr4A", "chr13", "chr14", "chr20", "chr15", "chr17", "chr18", "chr19", "chr21", "chr24", "chr26", "chr23", "chr28", "chr27", "chr22", "chr25", "chr29", "chrZ")
fst.labels <- c("2", "1", "3", "4", "1A", "5", "6","7","8","9","10", "11", "12", "4A", "13", "14", "20", "15", "17", "18", "19", "21", "24", "26", "23", "28", "27", "22", "25", "29", "Z")
fst.Z.label <- c("Z")

###
### Fst for P. atropygialis [PCA, N = 48] against P. cincta [PCC, N = 42]
###
fst1 <- read_delim("wg.pca_pcc.20kb_10kb_step.windowed.weir.fst", delim = "\t", col_names = c("CHROM", "START", "END", "VARIANTS", "WEIGHTED_FST", "MEAN_FST"), skip = 1)
fst1$window.id = paste(fst1$CHROM,fst1$START, sep=".")
fst1$CHROM=match(fst1$CHROM, order)
fst1=fst1[order(fst1$CHROM),]
## Evaluate the Fst landscape genome-wide, for the autosomes, and for chrZ
summary(fst1$WEIGHTED_FST)
fst1.auto <- subset(fst1, CHROM != "31")
summary(fst1.auto$WEIGHTED_FST)
quantile(fst1$WEIGHTED_FST, c(0.99), na.rm = TRUE)
fst1.Z <- subset(fst1, CHROM == "31")
summary(fst1.Z$WEIGHTED_FST)
quantile(fst1.Z$WEIGHTED_FST, c(0.99), na.rm = TRUE)

## Plot Fst in 20kb sliding windows
manhattan(fst1, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila atropygialis x P. cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(fst1.Z, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.Z.label, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila atropygialis x P. cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))


###
### Fst for P. atropygialis [PCA, N = 48] against P. hecki [PAH, N = 110]
###
fst2 <- read_delim("wg.pca_pah.20kb_10kb_step.windowed.weir.fst", delim = "\t", col_names = c("CHROM", "START", "END", "VARIANTS", "WEIGHTED_FST", "MEAN_FST"), skip = 1)
fst2$window.id = paste(fst2$CHROM,fst2$START, sep=".")
fst2$CHROM=match(fst2$CHROM, order)
fst2=fst2[order(fst2$CHROM),]
## Evaluate the Fst landscape genome-wide, for the autosomes, and for chrZ
summary(fst2$WEIGHTED_FST)
fst2.auto <- subset(fst2, CHROM != "31")
summary(fst2.auto$WEIGHTED_FST)
quantile(fst2.auto$WEIGHTED_FST, c(0.99), na.rm = TRUE)
fst2.Z <- subset(fst2, CHROM == "31")
summary(fst2.Z$WEIGHTED_FST)
quantile(fst2.Z$WEIGHTED_FST, c(0.99), na.rm = TRUE)

## Plot Fst in 20kb sliding windows
manhattan(fst2, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila a. hecki x P. c. atropygialis (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(fst2.Z, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.Z.label, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila a. hecki x P. c. atropygialis (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))


###
### Fst for P. cincta [PCC, N = 42] against P. hecki [PAH, N = 110]
###
fst3 <- read_delim("wg.pcc_pah.20kb_10kb_step.windowed.weir.fst", delim = "\t", col_names = c("CHROM", "START", "END", "VARIANTS", "WEIGHTED_FST", "MEAN_FST"), skip = 1)
fst3$window.id = paste(fst3$CHROM,fst3$START, sep=".")
fst3$CHROM=match(fst3$CHROM, order)
fst3=fst3[order(fst3$CHROM),]
## Evaluate the Fst landscape genome-wide, for the autosomes, and for chrZ
summary(fst3$WEIGHTED_FST)
fst3.auto <- subset(fst3, CHROM != "31")
summary(fst3.auto$WEIGHTED_FST)
quantile(fst3.auto$WEIGHTED_FST, c(0.99), na.rm = TRUE)
fst3.Z <- subset(fst3, CHROM == "31")
summary(fst3.Z$WEIGHTED_FST)
quantile(fst3.Z$WEIGHTED_FST, c(0.99), na.rm = TRUE)

## Plot Fst in 20kb sliding windows
manhattan(fst3, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila a. hecki x P. c. cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(fst3.Z, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.Z.label, bp="START", p="WEIGHTED_FST", logp=FALSE, ylab="Fst", xlab="Poephila a. hecki x P. c. cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))


######################
###Prepare PBS Data###
######################

## We are going to derive and plot the 'normalized population branch statistic' [PBSn1]
## See Equation #3 in Shpak et al. (2024) and Malaspinas et al. (2016)
## Note: Need to use Equations #1 and #2 from Shpak et al. (2024) too

# Find common positions across the three data frames
common_positions <- intersect(intersect(fst1$window.id, fst2$window.id), fst3$window.id)

fst1_filter <- fst1 |> filter(window.id %in% common_positions)
fst2_filter <- fst2 |> filter(window.id %in% common_positions)
fst3_filter <- fst3 |> filter(window.id %in% common_positions)

pbs <- fst1_filter[1:4]
pbs$window.id <- fst1_filter$window.id
#Equation #1: Distance function
pbs$Tac <- -log(1 - fst1_filter$WEIGHTED_FST)
pbs$Tah <- -log(1 - fst2_filter$WEIGHTED_FST)
pbs$Tch <- -log(1 - fst3_filter$WEIGHTED_FST)
#Equation #2: PBS
pbs$atro <- (pbs$Tah + pbs$Tac - pbs$Tch) / 2
pbs$cincta <- (pbs$Tch + pbs$Tac - pbs$Tah) / 2
pbs$hecki <- (pbs$Tah + pbs$Tch - pbs$Tac) / 2
#Equation #3: Normalized PBS [i.e., PBSn1]
pbs$n1a <- pbs$atro / (1 + pbs$atro + pbs$cincta + pbs$hecki)
pbs$n1c <- pbs$cincta / (1 + pbs$atro + pbs$cincta + pbs$hecki)
pbs$n1h <- pbs$hecki / (1 + pbs$atro + pbs$cincta + pbs$hecki)

## Evaluate which regions are more differentiated than others
#PBS
quantile(pbs$atro, c(0.999), na.rm = TRUE) #0.01% threshold: 0.4473701
pbs.atro.top99 <- pbs[pbs$atro >= 0.4473701,]
quantile(pbs$cincta, c(0.999), na.rm = TRUE) #0.01% threshold: 0.4405492
pbs.cincta.top99 <- pbs[pbs$cincta >= 0.4405492,]
#PBSn1
quantile(pbs$n1a, c(0.999), na.rm = TRUE) #0.01% threshold: 0.1436372
pbsn1.atro.top99 <- pbs[pbs$n1a >= 0.1436372,]
quantile(pbs$n1c, c(0.999), na.rm = TRUE) #0.01% threshold: 0.1528927
pbsn1.cincta.top99 <- pbs[pbs$n1c >= 0.1528927,]

## Are regions that have experienced greatest evolution correlated in each ssp?
cor.test(pbs$n1a, pbs$n1c, method=c("pearson")) #-0.1967188

## Find common positions across the three data frames
common_pbs_outliers <- intersect(pbs.atro.top99$window.id, pbs.cincta.top99$window.id)
common_pbsn1_outliers <- intersect(pbsn1.atro.top99$window.id, pbsn1.cincta.top99$window.id)

## Save output for evaluation
write.csv(pbsn1.atro.top99, file = "pbsn1.atro.top99.csv", row.names = FALSE)
write.csv(pbsn1.cincta.top99, file = "pbsn1.cincta.top99.csv", row.names = FALSE)

## Plot results
fst.labels <- c("2", "1", "3", "4", "1A", "5", "6","7","8","9","10", "11", "12", "4A", "13", "14", "20", "15", "17", "18", "19", "21", "24", "26", "23", "28", "27", "22", "25", "29", "Z")
#atropygialis
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="atro", logp=FALSE, ylab="PBS", xlab="Poephila atropygialis (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="n1a", logp=FALSE, ylab="Normalized PBS", xlab="Poephila atropygialis (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
#cincta
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="cincta", logp=FALSE, ylab="PBS", xlab="Poephila cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="n1c", logp=FALSE, ylab="Normalized PBS", xlab="Poephila cincta (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
#hecki
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="hecki", logp=FALSE, ylab="PBS", xlab="Poephila hecki (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))
manhattan(pbs, cex=0.35, cex.axis = 0.8, snp="window.id", chr="CHROM", chrlabs = fst.labels, bp="START", p="n1h", logp=FALSE, ylab="Normalized PBS", xlab="Poephila hecki (20kb - 10kb steps)", suggestiveline = F, genomewideline = F, ylim=c(0.0, 1.05))

plot(pbs.atro.top99$atro, pbs.atro.top99$n1a)

## Plot PBS in genomic windows with clusters of the most differentiated SNPs
pbs.subset <- pbs |> select(CHROM, window.id, START, END,  n1a, n1c)

## Chromosome 8: PTPRC
pbs.chr8 <- pbs.subset |> filter(CHROM == 9 & START >= 21550001 & END <= 22000001) |> pivot_longer(cols = c(n1a, n1c), names_to = "taxon", values_to = "PBS")
pbs.chr8 |> ggplot(aes(x = START, y = PBS, color = taxon)) +
  geom_line(size = 1.5, alpha = 0.75) +
  scale_color_manual(values=taxon_colors) +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  geom_segment(xend = 21796059, x = 21799285, y = 0.31, yend = 0.31, colour = "black") +
  geom_segment(x = 21626319, xend = 21688124, y = 0.25, yend = 0.25, colour = "black", linetype = "dashed", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 21724740, x = 21782596, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = 21805762, xend = 21817046, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 21870756, x = 21939306, y = 0.25, yend = 0.25, colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 21979488, x = 21994339, y = 0.25, yend = 0.25, colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  coord_cartesian(ylim = c(-0.05, 0.325)) +
  labs(x = "Chromosome 8", y = "PBS (normalized)") +
  theme(legend.position="none")

## Chromosome Z
pbs.chrZ <- pbs.subset |> filter(CHROM == 31) |> pivot_longer(cols = c(n1a, n1c), names_to = "taxon", values_to = "PBS")
  
pbs.chrZ |> filter(START >= 37500001 & END <= 39500001) |> ggplot(aes(x = START, y = PBS, color = taxon)) +
  geom_line(size = 1, alpha = 0.75) +
  scale_color_manual(values=taxon_colors) +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  geom_segment(xend = 38152507, x = 38168197, y = 0.31, yend = 0.31, colour = "black") +
  geom_segment(x = 38185813,  xend = 38211402, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 38247804,  x = 38255951, y = 0.25, yend = 0.25, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 38302682,  x = 38314521, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = 38314389,  xend = 38363799, y = 0.25, yend = 0.25, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = 38370056,  xend = 38374143, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = 39445857,  xend = 40160719, y = 0.25, yend = 0.25, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  coord_cartesian(ylim = c(-0.05, 0.325)) +
  labs(x = "Chromosome Z", y = "PBS (normalized)") +
  theme(legend.position="none")

pbs.chrZ |> filter(START >= 46500001 & END <= 48250001) |> ggplot(aes(x = START, y = PBS, color = taxon)) +
  geom_line(size = 1, alpha = 0.75) +
  scale_color_manual(values=taxon_colors) +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  geom_segment(xend = 46976005, x = 47074828, y = 0.31, yend = 0.31, colour = "black") +
  geom_segment(xend = 47492190, x = 47564810, y = 0.31, yend = 0.31, colour = "black") +
  geom_segment(xend = 48058877, x = 48073350, y = 0.31, yend = 0.31, colour = "black") +
  geom_segment(x = 46973557,  xend = 47016541, y = 0.3, yend = 0.3, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(xend = 46945287,  x = 46973519, y = 0.25, yend = 0.25, colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  coord_cartesian(ylim = c(-0.05, 0.325)) +
  labs(x = "Chromosome Z", y = "PBS (normalized)") +
  theme(legend.position="none")

pbs.chrZ |> filter(START >= 64500001 & END <= 66500001) |> ggplot(aes(x = START, y = PBS, color = taxon)) +
  geom_line(size = 1, alpha = 0.75) +
  scale_color_manual(values=taxon_colors) +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  geom_segment(xend = 65493363, x = 65568514, y = 0.31, yend = 0.31, colour = "black") +
  coord_cartesian(ylim = c(-0.05, 0.325)) +
  labs(x = "Chromosome Z", y = "PBS (normalized)") +
  theme(legend.position="none")


###############################
##IDENTIFYING CANDIDATE GENES##
###############################

## Read in gff for chromosome of interest - note that file has already been filtered to only include genes
setwd("/Users/danielhooper/Documents/Projects/LTF.Projects/reference/genome_assemblies_genome_gff")
gff8 <- read.table("GCF_003957565.2_bTaeGut1.4.pri_genomic.chr8.gff", sep="\t", header=F)
gffZ <- read.table("GCF_003957565.2_bTaeGut1.4.pri_genomic.chrZ.gff", sep="\t", header=F)

## Subset and clear up the gff - add names
colnames(gff8) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
colnames(gffZ) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

## Remove lncRNAs
gff8 <- gff8 |> filter(!grepl("lncRNA", attribute))
gffZ <- gffZ |> filter(!grepl("lncRNA", attribute))

## Arrange the gff
gff8 <- gff8 |> arrange(start, end)
gffZ <- gffZ |> arrange(start, end)

## Make a gene mid point variable
gff8 <- gff8 |> mutate(mid = start + (end-start)/2)
gffZ <- gffZ |> mutate(mid = start + (end-start)/2)
## Make a variable that can used to plot genes
gff8 <- gff8 |> mutate(plot = 1.2*(start/start))
gffZ <- gffZ |> mutate(plot = 1.2*(start/start))

## Identify the 10 most significant loci on a chromosome
hits <- pbs |> arrange(desc(n1a)) |> head(n = 10)

## Find the nearest genes to our highest hit
x <- hits$START[2]

## Find the set of genes within 150 kb of our top GWAS hit
gene_hits <- gff8 |> mutate(hit_dist = abs(mid - x)) |> arrange(hit_dist) |> filter(hit_dist < 100000)
gene_hits <- gffZ |> mutate(hit_dist = abs(mid - x)) |> arrange(hit_dist) |> filter(hit_dist < 100000)

## What are these genes?
gene_hits <- gene_hits |> dplyr::select(chr, start, end, attribute, hit_dist)

## Separate out the attribute column
gene_hits |> pull(attribute)
gene_hits |> pull(hit_dist)
