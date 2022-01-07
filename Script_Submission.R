#### Analysis of bacteria in feather-mosses phyllosphere ####
# Article title: Dominance of coniferous and broadleaved trees drives bacterial associations with boreal feather mosses
# Article authors: Juanita C. Rodríguez-Rodríguez, Yves Bergeron, Steven W. Kembel and Nicole J. Fenton
# Script written by: Juanita C. Rodríguez-Rodríguez

# 16S rRNA gene amplicon sequencing with universal primers 515F/926R
# 144 samples + 3 control (DNA extraction Kits) + 2 negative controls (PCR) + 1 positive control (PCR) = 149 samples sequenced
# Written on July 2020

# set file paths, # CHANGE ME to the directory containing the fastq files after unzipping.
path <- "/data/users/juanita/All_byPrimer/AllCyanobacteria/Juanita_Bryo 515-926"
list.files(path)

#Packages

### For DADA2 Analyses
#if (!requireNamespace("BiocManager", quietly = TRUE)) # https://benjjneb.github.io/dada2/dada-installation.html
#install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.11")
require(dada2) # https://github.com/benjjneb/dada2/issues/ If I have problems with dada2

### For Phyloseq and further statistical analyses
library(ggplot2); packageVersion("ggplot2")
library(picante); packageVersion("picante")
library(vegan)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")

library(devtools)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("seqTools") # Falta
library(seqtools)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

library(stats)
library(ape)
#library(ade4)

library(dplyr)
library(tidyr)

#library(adespatial) # For beta.div
library(nlme)
library(stringr) # For str_detect
library(multcomp) # For glht function to get comparisons from anova results
library(lsmeans) # To use lsmeans for ANOVA comparisons

#### DADA2 Analysis ####
# ---
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.", full.names = TRUE)) #filename of forward sequences
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.", full.names = TRUE)) #filename of reverse sequences

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plot quality profile of forward reads
plotQualityProfile(fnFs[1:4])
# Plot quality profile of reverse reads
plotQualityProfile(fnRs[1:4])

# need to trim 20 from start of both reads (gets rid of primer)
# Forward reads quality crashes around 270
# Reverse reads quality crashes around 220

# set filtered file folder path
filt_path <- file.path(path, "dada2-filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#### (1) TRIMMING AND FILTERING
# ---
# Filter reads at quality crashpoints identified earlier
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,20), truncLen=c(270,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=6, verbose=TRUE)

# Learn error rates
errF <- learnErrors(filtFs, multithread=6) # 107726250 total bases in 430905 reads from 17 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=6) # 100071600 total bases in 500358 reads from 21 samples will be used for learning the error rates.

# sanity check - visualize error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate the filtered fastq files (to analyze just unique sequences)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### (2) INFER SEQUENCE VARIANTS in each sample
# ---
dadaFs <- dada(derepFs, err=errF, pool="pseudo", multithread=6)
dadaRs <- dada(derepRs, err=errR, pool="pseudo", multithread=6)

# e.g. inspect results (dada-class: object describing DADA2 denoising results)
dadaFs[[1]] # 952 sequence variants were inferred from 8266 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaRs[[1]] # 775 sequence variants were inferred from 8221 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=15, returnRejects=FALSE, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### (3) CONSTRUCT SEQUENCE TABLE (is a higher-resolution analogue to the common OTU table)
# ---
seqtab <- makeSequenceTable(mergers)
sum(seqtab) # 1679342
dim(seqtab) # 149 30463 ASVs
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=6, verbose=TRUE) # Identified 14193 bimeras out of 23168 input sequences.
sum(seqtab.nochim) # 1423140
dim(seqtab.nochim) # 149 10870
sum(seqtab.nochim)/sum(seqtab) # 0.8474391

# summary - track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

sum(track[,1]) #To get the initial number of sequences (reads) produces by Illumina sequencing initially: 4297585

#### (4) TAXONOMICAL ASSIGMENT
# ---
# For 16S and universal primers is silva_nr_v132_train_set.fa.gz 
# I uses SILVA/UNITE for higher-level taxonomy (so for assignTaxonomy) and then assignSpecies with RefSeq+RDP for species-level.
taxa <- assignTaxonomy(seqtab.nochim, "/data/users/juanita/All_byPrimer/AllCyanobacteria/Juanita_Bryo 515-926/silva_nr_v132_train_set.fa.gz", multithread=6, tryRC=TRUE)

taxa.plus <- assignSpecies(taxa, "/data/users/juanita/All_byPrimer/AllCyanobacteria/Juanita_Bryo 515-926/RefSeq-RDP_dada2_assignment_species.fa.gz", allowMultiple = TRUE, tryRC = TRUE) #I should compare RefSeq with GTDB to see whichone is better

taxa.all <- as.matrix(data.frame(taxa, Genus_sp = paste(taxa.plus[,1], taxa.plus[,2], sep="_")))

taxa.print <- taxa.all # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
head(taxa.print, n=100)

dim(taxa) # 10870    6
dim(taxa.plus) # 10870    2
dim(taxa.all) # 10870    7
dim(seqtab.nochim) # 149 10870
dim(metadata) # 149   7


#### DATA Analysis with PHYLOSEQ ####
# ---
# Load metadata
head(sample.names)
metadata <- read.csv2("//data/users/juanita/All_byPrimer/AllCyanobacteria/Juanita_Bryo 515-926/Tmetadata.csv")
rownames(metadata) <- metadata$SampleID

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa.all))
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 10870 taxa and 149 samples ]
# sample_data() Sample Data:       [ 149 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 10870 taxa by 6 taxonomic ranks ]
sum(taxa_sums(ps)) #1423140
sum(sample_sums(ps)) #1423140

# Explore data
nsamples(ps) # 149
ntaxa(ps) # 10870 ASVs
sample_variables(ps)
rank_names(ps)
sample_data(ps)$Canopy
sample_data(ps)$Bry
subset_samples(ps, sample_data(ps)$Canopy == "BS")
subset_samples(ps, sample_data(ps)$Bry == "C")
subset_samples(ps, sample_data(ps)$sample_or_control == "control") # 5 samples
#number of seq per sample
summary(sample_sums(ps))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))]))
head(sort(sample_sums(ps),TRUE))
hist(sample_sums(ps))
#ASV richness
summary(estimate_richness(ps, measures = "Observed"))
#distribution of ASVs
hist(log10(taxa_sums(ps)))

### Working table filtering by  == Bacteria and excluding Chrloroplasts

ps.bacteria = subset_taxa(ps, Kingdom=="Bacteria")
ps.BnoC = subset_taxa(ps.bacteria, Order!="Chloroplast") #I took ps.BnoC for my community analysis (my metadata)
ps.BnoC
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 9409 taxa and 149 samples ]
# sample_data() Sample Data:       [ 149 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 9409 taxa by 7 taxonomic ranks ]
sum(taxa_sums(ps.BnoC)) #1177237 sequences
sum(sample_sums(ps.BnoC)) #1177237
class(ps.BnoC) # "phyloseq"

# CyanoBnoC
nsamples(ps.BnoC) # 149 samples
ntaxa(ps.BnoC) # 9409 ASVs / taxa

# Extract picante/vegan format objects
comm <- otu_table(ps.BnoC)
taxo <- data.frame(tax_table(ps.BnoC))
taxo$abund <- apply(comm, 2, sum)

# Combine the rest of metadata
metadata <- metadata[rownames(comm),]
head(metadata)

#Confirm if data matches
dim(comm)# 149 9409
dim(metadata) # 149   7
dim(taxo) # 9409    8

#To know if my names match
intersect(rownames(comm), rownames(metadata))
setdiff(rownames(comm), rownames(metadata))

# Abundant ASVs
taxo.print <- taxo
rownames(taxo.print) <- NULL
tail(taxo.print[order(taxo.print$abund),], n=20)

# Get abundances per taxon
aggregate(taxo$abund, by=list(taxo$Phylum), sum)
aggregate(taxo$abund, by=list(taxo$Class), sum)
aggregate(taxo$abund, by=list(taxo$Order), sum)
aggregate(taxo$abund, by=list(taxo$Family), sum)
aggregate(taxo$abund, by=list(taxo$Genus), sum)

# summary stats on communities and taxa
specnumber(comm)
hist(specnumber(comm))
hist(log10(apply(comm,1,sum)))
hist(log10(apply(comm,2,sum)))

#Exploring the data
#…how many samples was an ASV present in, that occurred more than 1 or 10 times?
dim(comm) # 149 9409
dim(comm[,apply(comm,2,sum)>1]) # 149 8725 ..occurred more than once
dim(comm[,apply(comm,2,sum)>10]) # 149 4200 ..occurred more thant 10 times
dim(comm[,apply(comm,2,sum)>100]) # 149 1439 ..occurred more thant 100 times
## Compare different cut offs

#…How many sequences per ASV are, that occurred in at least 2 samples and have more than 1 or 10 sequences? (gives the total number of sequences per ASV cutoff)
dim(comm[,apply(decostand(comm,method="pa"),2,sum)>1 & apply(comm,2,sum)>=100 ]) # 149 1437
dim(comm[,apply(decostand(comm,method="pa"),2,sum)>1 & apply(comm,2,sum)>1 ]) # 149 4533
dim(comm[,apply(decostand(comm,method="pa"),2,sum)>1 & apply(comm,2,sum)>=10 ]) # 149 3721 (3721 ASVs had atleast 11 sequences). I select this one: Ocurred in at least 2 samples and have at least 10 sequences. 
#For the analysis I select samples that are present in at least 2 ASVs (more than once) and also the ASVs that occurred in at least 2 samples and have at least 10 sequences.
#...but to start, follow the steps:

##STEPS: 
#1) Get rid of controls (I get comm.sub)
#2) Get rid of the rarests AVSs (I used: Numer of ASVs that ocurred in at least 2 samples and have more than 10 sequences) (I get comm.sub2)
#3) Going back to the same object, get rid of the ASVs that are 0
#4) Use this comm.sub2 for DESeq2 analysis
#5) Rarefy (I get comm.sub.rare)
#6) Do the rest of community alpha and beta analysis


### STEP 1) Checking for controls:
#2. Remove outliers - NMDS
#relative abundance
ps.ra <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ordinate
ps.ra.nmds <- ordinate(ps, method = "NMDS", k = 2, try = 100, distance = "bray") 
#Stress: Run 100 stress 9.459718e-05 
#... Procrustes: rmse 0.05221173  max resid 0.3167265 
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#87: stress < smin
#12: scale factor of the gradient < sfgrmin
plot_ordination(ps.ra, ps.ra.nmds, color = "Canopy", shape = "sample_or_control") +
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = SampleID), check_overlap = FALSE, size = 5) +
  geom_point(size = 1) + scale_shape_manual(values = c(19, 1))

#Remove outliers 1 (controls)
out.controls <-  c("PCR-CTRL-neg-Juanita","Control-Kit-1-Juanita","Control-Kit-2-Juanita","Control-Kit-3-Juanita","PCR-CTRL-Pos-Juanita")
ps.ra.noctl <- prune_samples(!sample_data(ps.ra)$SampleID %in% out.controls, ps.ra)
ps.ra.noctl <- prune_taxa(taxa_sums(ps.ra.noctl)>0,ps.ra.noctl)
ps1.ra <- transform_sample_counts(ps.ra.noctl, function(otu) otu/sum(otu)) 

ps.ra.nmds1 <- ordinate(ps1.ra, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps1.ra, ps.ra.nmds1, color = "Canopy", shape = "sample_or_control") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = SampleID), check_overlap = FALSE, size = 3, nudge_y = -0.1) + 
  geom_point(size = 5) + scale_shape_manual(values = c(19, 1))

apply(comm, 1, sum) # To know how many ASVs for each control are
# Positive control (20199), PCR-CTRL-neg-Juanita (10), Control-Kit-1-Juanita (62), Control-Kit-2-Juanita (19), Control-Kit-3-Juanita (129)

# Eliminate Positive and negative controls
comm.sub <- comm[-which(rownames(comm)=="PCR-CTRL-Pos-Juanita"),] #Here I eliminate Pos. control
comm.sub <- comm.sub[-which(rownames(comm.sub)=="Control-Kit-1-Juanita"),]
comm.sub <- comm.sub[-which(rownames(comm.sub)=="Control-Kit-2-Juanita"),]
comm.sub <- comm.sub[-which(rownames(comm.sub)=="Control-Kit-3-Juanita"),]
comm.sub <- comm.sub[-which(rownames(comm.sub)=="PCR-CTRL-neg-Juanita"),]

rownames(comm.sub)
metadata.sub <- metadata[rownames(comm.sub),] # match the data set with the others

#To know if my names match
intersect(rownames(comm.sub), rownames(metadata.sub))
setdiff(rownames(comm.sub), rownames(metadata.sub))

#Check both
setdiff(rownames(comm), rownames(comm.sub)) # I get the five controls
dim(comm) # 149 9409
dim(comm.sub) # 144 9409 (without controls)
dim(metadata.sub) # 144   7

### STEP 2) Get rid of the rarests AVSs (I used: Number of ASVs that occurred in at least 2 samples and have more than 10 sequences) (I get comm.sub2)
comm.sub2 <- comm.sub[,apply(decostand(comm.sub,method="pa"),2,sum)>1 & apply(comm.sub,2,sum)>=10 ] # >1 sequence per sample in, and >=10 sequences per ASVs
dim(comm.sub2) #144 3714 : There are 3714 sequences that occurred in at least 2 samples and have at least 10 sequences per ASV (total number of sequences per ASV cutoff)

### STEP 3) Going back to the same object and:
# a) Take the samples have at least 1000 sequences left:
comm.sub2 <- comm.sub2[apply(comm.sub2,1,sum)>=1000,] #To see how many samples have at least 1000 sequences left (this after get rid of the rearest ASVs)
# b) Get rid of the ASVs that have 0 sequences
comm.sub2 <- comm.sub2[,apply(comm.sub2,2,sum)>0]
taxo.sub <- taxo[colnames(comm.sub2),] #make them match with the names
taxo.sub$abund <- apply(comm.sub2,2,sum)

### STEP 4) Use this comm.sub2 for DESeq2 analysis
### Analysis DESeq2 (non-rarefied data)
comm.sub.de2 <- as.data.frame(comm.sub2 + 1) # get rid of zeroes by adding 1 to all abundances
de2dat <- DESeqDataSetFromMatrix(countData=as.matrix(t(comm.sub.de2)), colData=metadata.sub, design = ~ Canopy)
de2dat.de2 <- DESeq(de2dat)
#Get differential expression results
res <- results(de2dat.de2)
summary(res)

# plot LFCs
plotMA(res)

# combine results with taxonomy to ID interesting taxa
resOrdered <- res[order(res$padj),]
resOrdered <- resOrdered[which(resOrdered$padj<=0.05),]
resOrdered.taxo <- taxo.sub[rownames(resOrdered),]
resOrdered <- cbind(resOrdered, resOrdered.taxo)
# summarize - taxonomy of all OTUs with adjusted p-value <= 0.05
resOrdered.2 <- resOrdered
rownames(resOrdered.2) <- NULL
resOrdered.2 <- as.data.frame(resOrdered.2)

dim(resOrdered.2) # 1959   14
resOrdered[1,]

# example plots for some OTUs with greatest difference between Canopy (for the first line)
plotCounts(de2dat.de2 , gene="ACGGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCCGCAGGTGGCGAAGTAAGTCTGCTGTTAAAGCGTCTAGCTCAACTAGATAAGAGCAGTGGAAACTACTTACGCTAGAGTGCGTTCGGGGCAGAGGGAATTCCTGGTGTAGCGGTGAAATGCGTAGAGATCAGGAAGAACACCAGTGGCGAAGGCGCTCTGCTAGGCCGCAACTGACACTGAGGGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTAGCCGTAAACGATGGATACTAGGCGTGGCTTGTATCGACCCGAGCCGTGCCGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTACGCACGCAAGTGTG", intgroup="Canopy")

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x = tapply(resOrdered.2$log2FoldChange, resOrdered.2$Phylum, function(x) max(x))
x = sort(x, TRUE)
resOrdered.2$Phylum = factor(as.character(resOrdered.2$Phylum), levels=names(x))

# Family order
x = tapply(resOrdered.2$log2FoldChange, resOrdered.2$Family, function(x) max(x))
x = sort(x, TRUE)
resOrdered.2$Family = factor(as.character(resOrdered.2$Family), levels=names(x))

# !!! DESeq2 Graph of resOrdered.2 with Family
Fig_deseq <- ggplot(resOrdered.2, aes(x=Family, y=log2FoldChange, color=Phylum)) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  labs(title="Phyllosphere diversity of feather mosses under each canopy dominance",x="Bacterial Family", y = "log2FoldChange")+ 
  scale_color_manual(values=c('chartreuse3','gold','turquoise','coral','gold3','violet','palegreen','palevioletred','lightpink','turquoise4','mediumorchid3'))

#---#

### STEP 5) Rarefy

##EXTRA: Other checkings at different points:
# what's smallest number of sequences per sample (with at least...(No.))
min(apply(comm.sub2[apply(comm.sub2,1,sum)>1000,], 1, sum)) # 1920
min(apply(comm.sub2[apply(comm.sub2,1,sum)>2000,], 1, sum)) # 2029
min(apply(comm.sub2[apply(comm.sub2,1,sum)>2500,], 1, sum)) # 2720
min(apply(comm.sub2[apply(comm.sub2,1,sum)>3000,], 1, sum)) # 3095
min(apply(comm.sub2[apply(comm.sub2,1,sum)>4000,], 1, sum)) # 4016
#To check how many samples there are left:
dim(comm.sub2[apply(comm.sub2,1,sum)>1000,]) # 144 3714
dim(comm.sub2[apply(comm.sub2,1,sum)>2000,]) # 142 3714
dim(comm.sub2[apply(comm.sub2,1,sum)>2500,]) # 136 3714
dim(comm.sub2[apply(comm.sub2,1,sum)>3000,]) # 129 3714
dim(comm.sub2[apply(comm.sub2,1,sum)>4000,]) # 118 3714

#Now check which samples I loose comparing two cut offs:
setdiff(rownames(comm.sub2[apply(comm.sub2,1,sum)>1920,]), rownames(comm.sub2[apply(comm.sub2,1,sum)>3000,]) ) # "BE1-CD1" "BE1-CD5" "BE1-CD6" "BE1-SD5" "BE2-CD5" "BE2-CD6" "BE2-SD3" "BE2-SD5" "BE2-SD6" "BP2-SD5" "BP2-SD6" "CE1-CD4" "CE2-SD3" "CP1-SD4"
setdiff(rownames(comm.sub2[apply(comm.sub2,1,sum)>1920,]), rownames(comm.sub2[apply(comm.sub2,1,sum)>2000,]) ) # "BP2-SD5"
setdiff(rownames(comm.sub2[apply(comm.sub2,1,sum)>1920,]), rownames(comm.sub2[apply(comm.sub2,1,sum)>4000,]) ) # "BP2-SD5"
setdiff(rownames(comm.sub2[apply(comm.sub2,1,sum)>1920,]), rownames(comm.sub2[apply(comm.sub2,1,sum)>5000,]) ) # "BP2-SD5"

setdiff(rownames(comm.sub2[apply(comm.sub2,1,sum)>1920,]), rownames(comm.sub2[apply(comm.sub2,1,sum)>2500,]) ) # "BE1-CD1" "BE1-CD5" "BE1-CD6" "BE2-CD5" "BE2-SD5" "BP2-SD5" "BP2-SD6"

### Subset and Rarefy
### Continue with RAREFACTION analysis
# what's smallest number of sequences per sample (with at least 1000)?
min(apply(comm.sub2[apply(comm.sub2,1,sum)>1000,], 1, sum)) # 1920

rarecurve(comm.sub2, step=100, labe=FALSE)
abline(v=(min(rowSums(comm.sub2)))) # and adding a vertical line at the fewest seqs in any sample (correspond to the number obtained before = 1920)
dim(comm.sub2[apply(comm.sub2,1,sum)>1920,]) # 143 3714 ##To check how many samples there are left
comm.sub.rare <- rrarefy(comm.sub2, sample=1920) 
#This are taking samples with minimum 1920 sequences
comm.sub.rare <- comm.sub.rare[,apply(comm.sub.rare,2,sum)>0]
taxo.sub.rare <- taxo.sub[colnames(comm.sub.rare),] 
taxo.sub.rare$abund <- apply(comm.sub.rare,2,sum)
metadata.sub2 <- metadata.sub[rownames(comm.sub.rare),] # match the data set with the others
metadata.sub2$Canopy <- factor(metadata.sub2$Canopy)

dim(comm.sub.rare) #144 3694
sum(taxa_sums(comm.sub.rare)) #276480 (no. of sequences) 
sum(sample_sums(comm.sub.rare)) #276480 
sample_sums(comm.sub.rare) # To verify that I get the same number of sequences in each sample (1920), which is the cutoff of the rarefaction

### PERMANOVA
# !!! PERMANOVA: Test for community composition differences among Canopy
set.seed(25)
perm.comm.sub.rare <- adonis2(decostand(comm.sub.rare, method = "hellinger") ~ Canopy * Bry, data=metadata.sub2, permutations = how(blocks = metadata.sub2$Site, nperm = 9999), method = "bray")
perm.comm.sub.rare

### STEP 6) Do the rest of community alpha and beta analysis
#### STATISTICAL Analyses ####
# ---
# Analyses for all data, with rarefied to 1920 samples
# comm.sub.rare - Community data (Sequence table with abundances of ASVs in the different samples)
# taxo.sub.rare - OTU taxonomy data (Taxonomic IDs for the ASVs)
# metadata.sub2 - Metadata (All related information for each sample)
# ---

# Function taxocomm from https://github.com/skembel/seqtools/blob/master/R/taxocomm.R
taxocomm <- function(comm, taxo, rank) {
  comm.taxo <- aggregate(t(comm), by=list(taxonrank=taxo[,rank]), sum)
  rownames(comm.taxo) <- comm.taxo[,1]
  return(t(comm.taxo[,-1]))
}

#Seleccion to add to the ordiplot the arrows with "Phylum"
comm.phylum <- taxocomm(comm.sub.rare, taxo.sub.rare, "Phylum")
comm.phylum.rel <- decostand(comm.phylum, method="hellinger")
apply(comm.phylum.rel[,apply(comm.phylum.rel,2,mean)>0.01], 2, mean)

### Ordination
## !!! NMDS!
set.seed(27)
comm.sub.rare.mds <- metaMDS(comm.sub.rare) #distance = "bray" is the default
stressplot(comm.sub.rare.mds)
#Stress: 0.1121544 

ordi.nmds <- ordiplot(comm.sub.rare.mds, display="sites", type="points", cex = 0.2)
points(ordi.nmds$sites[metadata.sub2$Canopy == 'BS' & metadata.sub2$Bry == 'C',], pch=24, col="darkgreen", bg="chartreuse3", cex = 0.7)
points(ordi.nmds$sites[metadata.sub2$Canopy == 'TA' & metadata.sub2$Bry == 'C',], pch=21, col="darkgreen", bg="chartreuse3", cex = 0.7)
points(ordi.nmds$sites[metadata.sub2$Canopy == 'BS' & metadata.sub2$Bry == 'S',], pch=24, col="orangered4", bg="orangered1", cex = 0.7)
points(ordi.nmds$sites[metadata.sub2$Canopy == 'TA' & metadata.sub2$Bry == 'S',], pch=21, col="orangered4", bg="orangered1", cex = 0.7)
ordiellipse(comm.sub.rare.mds, metadata.sub2$Bry, label=TRUE, kind = c("sd"), col = "chartreuse3", cex=1.2, font=2, lty=2, lwd=3)
ordiellipse(comm.sub.rare.mds, metadata.sub2$Canopy, label=TRUE, kind = c("sd"), col = "dodgerblue4", cex=1.2, font=3, lty=1, lwd=3)
plot(envfit(ordi.nmds, comm.phylum.rel), p.max=0.05, col="black", cex=1)

envf <- envfit(comm.sub.rare.mds, metadata.sub2[,c("Canopy", "Bry", "Site","Block")])
scores.envfit

# summary stats on communities and taxa
specnumber(comm.sub.rare)
specnumber(comm.sub.rare, metadata.sub2$Canopy) # BS 2436, TA 3340 Gamma diversity
hist(specnumber(comm.sub.rare))

# Diversity
boxplot(diversity(comm.sub.rare) ~ paste(metadata.sub2$Canopy, metadata.sub2$Bry))
boxplot(diversity(comm.sub.rare) ~ paste(metadata.sub2$Canopy))

div.model_bryyy <- lme(diversity(comm.sub.rare, index = "shannon") ~ Canopy, random = ~1|Site/Block, data = metadata.sub2)
summary(div.model_bryyy)
anova(div.model_bryyy)
plot(div.model_bryyy)

div.model <- lme(diversity(comm.sub.rare, index = "shannon") ~ Canopy*Bry, random = ~1|Site/Block, data = metadata.sub2)
summary(div.model)
anova(div.model) # !!!
plot(div.model)
hist(diversity(comm.sub.rare, index = "shannon"))
length(diversity(comm.sub.rare, index = "shannon"))

hist(diversity(comm.sub.rare.phylum, index = "shannon"))

summary(glht(div.model, linfct=mcp(Canopy="Tukey")), test = adjusted(type = "bonferroni"))
lsmeans(div.model, pairwise~Canopy, adjust="tukey") # !!!
lsmeans(div.model, pairwise~Canopy*Bry, adjust="tukey") #I get all the contrasts

#Diversity
ggplot(metadata.sub2, aes(x = paste(metadata.sub2$Canopy, metadata.sub2$Bry), y = diversity(comm.sub.rare), fill=metadata.sub2$Canopy)) + geom_violin()

##Diversity grided by Canopy
Canopy_names <- c("BS"="Black spruce", "TA"="Trembling aspen")
Bry_names <- c("C"="P. crista-castrensis", "S"="P. schreberi")

#!!! Specnumber by Bry and Canopy (in boxplot)
div.boxplot <- ggplot(metadata.sub2, aes(x = metadata.sub2$Bry, y = diversity(comm.sub.rare, index = "shannon"), fill=metadata.sub2$Bry)) +
  geom_boxplot(fill="white") +
  labs(title="Variation of feather-mosses phyllosphere between Sites",x="Sites", y = "Specnumber") +
  facet_grid(cols=vars(metadata.sub2$Canopy), labeller = as_labeller(Canopy_names))

div.boxplot + scale_fill_manual(values=alpha(c("chartreuse3" ,"orangered1", "chartreuse3" ,"orangered1"), 0.7), name = "Feather-mosses") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + theme(strip.background = element_rect(colour = "black", fill = "white"))

#!!! Specnumber by Canopy (in boxplot)
div.boxplot2 <- ggplot(metadata.sub2, aes(x = metadata.sub2$Canopy, y = diversity(comm.sub.rare, index = "shannon"), fill=metadata.sub2$Canopy)) +
  geom_boxplot(fill="white") +
  labs(title="Variation of feather-mosses phyllosphere between Sites",x="Sites", y = "Shannon")

div.boxplot2 + scale_fill_manual(values=alpha(c("chartreuse3" ,"orangered1", "chartreuse3" ,"orangered1"), 0.7), name = "Feather-mosses") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + theme(strip.background = element_rect(colour = "black", fill = "white"))


#!!! Diversity plot with Shannon index
div.plot <- ggplot(metadata.sub2, aes(x = metadata.sub2$Bry, y = diversity(comm.sub.rare, index = "shannon"), fill=metadata.sub2$Bry)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Phyllosphere diversity of feather mosses under each canopy dominance",x="Feather-mosses", y = "Shannon Diversity")+
  facet_grid(cols=vars(metadata.sub2$Canopy), labeller = as_labeller(Canopy_names))
# Use custom color palettes
div.plot + scale_fill_manual(values=alpha(c("chartreuse3" ,"orangered1", "chartreuse3" ,"orangered1"), 0.7), name = "Feather-mosses") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + theme(strip.background = element_rect(colour = "black", fill = "white"))

#Specnumber (finds the No. of spp, With MARGIN = 2, it finds frequencies of species)
div.spec <- ggplot(metadata.sub2, aes(x = metadata.sub2$Bry, y = specnumber(comm.sub.rare), fill=metadata.sub2$Bry)) + geom_violin() + geom_boxplot(width=0.1, fill="white")+
  labs(title="Phyllosphere diversity of feather mosses under each canopy dominance",x="Feather-mosses", y = "Diversity")+
  facet_grid(cols=vars(metadata.sub2$Canopy), labeller = as_labeller(Canopy_names))
# Use custom color palettes
div.spec + scale_fill_manual(values=alpha(c("chartreuse3" ,"orangered1", "chartreuse3" ,"orangered1"), 0.7), name = "Feather-mosses") + theme(panel.background = element_rect(fill = "white", colour = "grey50")) + theme(strip.background = element_rect(colour = "black", fill = "white"))

#!!! Specnumber by Site (in boxplot)
spec.sites2 <- ggplot(metadata.sub2, aes(factor(metadata.sub$Site, levels=c("A","B","C")), y=specnumber(comm.sub.rare), fill=metadata.sub2$Site)) +
  geom_boxplot(fill="white")+
  labs(title="Variation of feather-mosses phyllosphere between Sites",x="Sites", y = "Specnumber") + facet_grid(cols=vars(metadata.sub2$Canopy))

spec.sites2 + scale_fill_manual(values=alpha(c("slategray1" ,"slategray3", "slategray4", "slategray1" ,"slategray3", "slategray4"), 0.7), name = "Feather-mosses") + theme(panel.background = element_rect(fill = "white", colour = "grey50")) + theme(strip.background = element_rect(colour = "black", fill = "white"))

# Alpha diversity
plot_richness(comm.sub.rare, x=metadata.sub2$Canopy, measures=c("Shannon", "Simpson"), color=metadata.sub2$Bry) #???

##Beta diversity (post test)
spe.dist<-vegdist(comm.sub.rare, method="bray")
#spe.dist.perm <- adonis(spe.dist ~metadata.sub2$Canopy, data=metadata.sub2,perm=9999)

dispersion <- betadisper(spe.dist, group=metadata.sub2$Canopy,type = c("median","centroid"), bias.adjust = TRUE)
dispersion # !!! To get the contrast (average distance to median)
permutest(dispersion) # !!! TO get the P value
summary(dispersion)
pair.adonis <- pairwise.adonis(spe.dist,metadata.sub2$Canopy,perm=9999) #Function installed from https://github.com/pmartinezarbizu/pairwiseAdonis

dispersion.Bry <- betadisper(spe.dist, group=metadata.sub2$Bry,type = c("median","centroid"), bias.adjust = TRUE)
dispersion.Bry
permutest(dispersion.Bry)

### Who are the most abundant ASVs and taxa?
# list most abundant ASVs
taxo.sub.rare2 <- taxo.sub.rare
rownames(taxo.sub.rare2) <- NULL #Because otherwise I get the whole sequence names.
head(taxo.sub.rare2)
head(taxo.sub.rare2, n=100)

taxo.sub.rare2$abund <- apply(comm.sub.rare, 2, sum)
taxo.sub.rare2$relabund <- taxo.sub.rare2$abund / sum(taxo.sub.rare2$abund)
tail(taxo.sub.rare2[order(taxo.sub.rare2$relabund),], n=30)

## Aggregate community at different taxonomic levels
comm.sub.rare.phylum <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Phylum")
comm.sub.rare.class <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Class")
comm.sub.rare.order <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Order")
comm.sub.rare.family <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Family")
comm.sub.rare.genus <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Genus")
comm.sub.rare.species <- taxocomm(comm.sub.rare, taxo.sub.rare2, "Genus_sp")

# list most abundant taxa at different taxonomic ranks (No. of seq per category)
apply(comm.sub.rare.phylum,2,sum)[order(apply(comm.sub.rare.phylum,2,sum))]
apply(comm.sub.rare.class,2,sum)[order(apply(comm.sub.rare.class,2,sum))]
apply(comm.sub.rare.order,2,sum)[order(apply(comm.sub.rare.order,2,sum))]
apply(comm.sub.rare.family,2,sum)[order(apply(comm.sub.rare.family,2,sum))]
apply(comm.sub.rare.genus,2,sum)[order(apply(comm.sub.rare.genus,2,sum))]
apply(comm.sub.rare.species,2,sum)[order(apply(comm.sub.rare.species,2,sum))]

# Get percentage of all sequences for Phylum in both canopies...
apply(comm.sub.rare.phylum,2,sum)[order(apply(comm.sub.rare.phylum,2,sum))]/sum(comm.sub.rare.phylum)

### Summarize community taxonomic structure plots
# Phylum - all samples
comm.sub.rare.phylum.long <- matrix2sample(comm.sub.rare.phylum)
ggplot(comm.sub.rare.phylum.long, aes(x=plot, y=abund, fill=id)) + geom_bar(stat="identity")
# Family - all samples
comm.sub.rare.family.long <- matrix2sample(comm.sub.rare.family)
ggplot(comm.sub.rare.family.long, aes(x=plot, y=abund, fill=id)) + geom_bar(stat="identity")
# Class - all samples
comm.sub.rare.class.long <- matrix2sample(comm.sub.rare.class)
ggplot(comm.sub.rare.class.long, aes(x=plot, y=abund, fill=id)) + geom_bar(stat="identity")

# Phylum - organize per Canopy
comm.sub.rare.phylum.bytrt <- aggregate(comm.sub.rare.phylum, by=list(type=metadata.sub2$Canopy), mean)
rownames(comm.sub.rare.phylum.bytrt) <- comm.sub.rare.phylum.bytrt[,1]
comm.sub.rare.phylum.bytrt <- comm.sub.rare.phylum.bytrt[,-1]
comm.sub.rare.phylum.bytrt <-  decostand(comm.sub.rare.phylum.bytrt, method="hellinger")
comm.sub.rare.phylum.bytrt.long <- matrix2sample(comm.sub.rare.phylum.bytrt)

# Class - organize per Canopy
comm.sub.rare.class.bytrt <- aggregate(comm.sub.rare.class, by=list(type=metadata.sub2$Canopy), mean)
rownames(comm.sub.rare.class.bytrt) <- comm.sub.rare.class.bytrt[,1]
comm.sub.rare.class.bytrt <- comm.sub.rare.class.bytrt[,-1]
comm.sub.rare.class.bytrt <-  decostand(comm.sub.rare.class.bytrt, method="hellinger")
comm.sub.rare.class.bytrt.long <- matrix2sample(comm.sub.rare.class.bytrt)

# Sort by mean of Phylum value to make plot look better
comm.sub.rare.phylum.bytrt.long$id <- factor(comm.sub.rare.phylum.bytrt.long$id, levels=colnames(comm.sub.rare.phylum.bytrt)[order(apply(comm.sub.rare.phylum.bytrt,2,mean), decreasing = TRUE)])

# Sort by mean of Class value to make plot look better
comm.sub.rare.class.bytrt.long$id <- factor(comm.sub.rare.class.bytrt.long$id, levels=colnames(comm.sub.rare.class.bytrt)[order(apply(comm.sub.rare.class.bytrt,2,mean), decreasing = TRUE)])

#!!! Stacked lineplot of Phylum by Canopy
Stacked_Canopy <- ggplot(comm.sub.rare.phylum.bytrt.long, aes(x=plot, y=abund, group=id, color=id)) + geom_point() + geom_line(aes(color=id), size=1.5)+
  labs(title="Phyla of feather-mosses phylloshere", x="Forest type", y="Relative abundance")+
  theme_bw()+
  scale_color_manual(values=c('turquoise','coral','gold','chartreuse3','violet','gold3','lightpink','palegreen','turquoise4','mediumorchid3','olivedrab1','thistle','palevioletred','slateblue1','darkorange','snow4'))

### !!! Test ANOVA for differences in taxon abundances between Canopies for Phyllum Cyanobacteria
##This is the original loop (witouth calculatiing p.adj.)
# Phylum from log10 data
for (taxon in colnames(comm.sub.rare.phylum)) {
  print(taxon)
  print(anov <- anova(lm(log10(comm.sub.rare.phylum[,taxon]+1) ~ metadata.sub2$Canopy)))
}

## Get p.adjust of thep values for each Phyla
P.Table<-data.frame(NULL)

for (i in 1:dim(comm.sub.rare.phylum)[2]) {
  
  p.values <- data.frame(Taxon= colnames(comm.sub.rare.phylum)[i],
                         Pvalue = round(anova(lm(log10(comm.sub.rare.phylum[,i]+1) ~ metadata.sub2$Canopy))$`Pr(>F)` [1], 5))
  
  # adds the row for each Taxon (p and p adjusted)
  P.Table<-rbind(P.Table,p.values )  }

## claculating P.adj
P.Table$P.adj <- p.adjust(P.Table$Pvalue , "BH")
P.Table

### Same graph Stacked lineplot but for Bry and the ANOVA test
## Prepare comm.sub.rare.phylum.bytrt.forBry.long
# Phylum - organize per Bry
comm.sub.rare.phylum.bytrt.forBry <- aggregate(comm.sub.rare.phylum, by=list(type=metadata.sub2$Bry), mean)
rownames(comm.sub.rare.phylum.bytrt.forBry) <- comm.sub.rare.phylum.bytrt.forBry[,1]
comm.sub.rare.phylum.bytrt.forBry <- comm.sub.rare.phylum.bytrt.forBry[,-1]
comm.sub.rare.phylum.bytrt.forBry <-  decostand(comm.sub.rare.phylum.bytrt.forBry, method="hellinger")
comm.sub.rare.phylum.bytrt.forBry.long <- matrix2sample(comm.sub.rare.phylum.bytrt.forBry)

# Sort by mean of Phylum value to make plot look better
comm.sub.rare.phylum.bytrt.forBry.long$id <- factor(comm.sub.rare.phylum.bytrt.forBry.long$id, levels=colnames(comm.sub.rare.phylum.bytrt.forBry)[order(apply(comm.sub.rare.phylum.bytrt.forBry,2,mean), decreasing = TRUE)])

#!!! Stacked lineplot of Phylum by Bryophyte
Stacked_Bry <- ggplot(comm.sub.rare.phylum.bytrt.forBry.long, aes(x=plot, y=abund, group=id, color=id)) + geom_point() + geom_line(aes(color=id), size=1.5)+
  labs(title="Bacterial phylla associated with feather-mosses", x="Feather-moss species", y="Relative abundance")+
  theme_bw()+
  scale_color_manual(values=c('turquoise','coral','gold','chartreuse3','violet','gold3','lightpink','palegreen','turquoise4','mediumorchid3','olivedrab1','thistle','palevioletred','slateblue1','darkorange','snow4')) 
Stacked_Bry

### !!! ANOVA for differences in taxon abundances between Byr for Phyllum Cyanobacteria
##This is the original loop (witouth calculatiing p.adj.
# Phylum from log10 data
for (taxon in colnames(comm.sub.rare.phylum)) {
  print(taxon)
  print(anova(lm(log10(comm.sub.rare.phylum[,taxon]+1) ~ metadata.sub3$Bry)))
}

## Get p.adjust of thep values for each Phyla for Bry
P.Table.bry<-data.frame(NULL)

for (i in 1:dim(comm.sub.rare.phylum)[2]) {
  
  p.values <- data.frame(Taxon= colnames(comm.sub.rare.phylum)[i],
                         Pvalue = round(anova(lm(log10(comm.sub.rare.phylum[,i]+1) ~ metadata.sub3$Bry))$`Pr(>F)` [1], 5))
  
  # adds the row for each Taxon (p and p adjusted)
  P.Table.bry<-rbind(P.Table.bry,p.values )  }

## claculating P.adj
P.Table.bry$P.adj <- p.adjust(P.Table.bry$Pvalue , "BH")
P.Table.bry

#!!! PERMANOVA: To test differences in relative abundance of Phyla (total abundance of the group)
set.seed(35)
perm.comm.phylum.rel2 <- adonis2(comm.phylum.rel ~ Canopy * Bry, data=metadata.sub2, permutations = how(blocks = metadata.sub2$Site, nperm = 9999), method = "euclidean") # I get sig. Canopy and Bry (0.0001) but not the interaction Canopy:Bry
perm.comm.phylum.rel2

#!!! ANOVA general for phyllum!!!
div.model.all <- lme(diversity(comm.phylum.rel) ~ Canopy*Bry, random = ~1|Site/Block, data = metadata.sub2)
anova(div.model.all)