library(rio)
library(dplyr)
library(tidyverse)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

col_agegroup <- c("#86B875","#4CB9CC","#CD99D8")
col_gender <- c("#E15566","#4981BF")

meta <- import(snakemake@input[["meta"]])
print(head(meta))

meta <- meta %>%
	filter(Cohort == "Johansson et al.") %>% dplyr::rename(Basename=barcode) %>% 
	select(-`Cell Type`, -`Age Group`) %>% 
	mutate(AgeGroup = ifelse(Age < 40, "<40", 
				ifelse(Age > 60, ">60", "40-60")))
print(head(meta))
print(table(meta$AgeGroup, meta$Gender))
export(meta, snakemake@output[["metadata"]])

rgSet <- read.metharray.exp(file.path("/home/users/nus/e0859928/scratch/Snakemake/pbmc-aging-snakemake/data/dna_methylation"), targets=meta)
sampleNames(rgSet) <- meta$Basename
annotation(rgSet)["array"] = "IlluminaHumanMethylation450k"
annotation(rgSet)["annotation"] = "ilmn12.hg19"

qcReport(rgSet, sampGroups=meta$AgeGroup, 
         pdf=snakemake@output[["qcplot"]])
saveRDS(rgSet, snakemake@output[["rgset"]])

## Normalization
mSetNoob <- preprocessNoob(rgSet, verbose = TRUE) 
print("here")
saveRDS(mSetNoob, snakemake@output[["mset"]])
print("here")

## PC1 vs. PC2
pdf(snakemake@output[["pcaplot"]], width = 12, height = 6)
par(mfrow=c(1,2))
plotMDS(getM(mSetNoob), top=1000, gene.selection="common", 
        col=col_agegroup[factor(meta$AgeGroup, levels = c("<40","40-60",">60"))])
legend("topright", legend=c("<40","40-60",">60"), text.col=col_agegroup,
       bg="white", cex=0.7)

plotMDS(getM(mSetNoob), top=1000, gene.selection="common",  
        col=col_gender[factor(meta$Gender, levels = c("female", "male"))])
legend("top", legend=c("female", "male"), text.col=col_gender,
       bg="white", cex=0.7)
dev.off()

