library(rio)
library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)
library(ensembldb)

count <- import(snakemake@input[["mtx"]])
metadata <- import(snakemake@input[["metadata"]])
colnames(metadata) <- c("cell_ID", "donor_id","disease")
metadata <- metadata %>% mutate(disease = ifelse(disease == "CT", "normal", "supercentenarian"))
count <- count %>% column_to_rownames("V1")

print(count[1:5,1:5])
print(unique(metadata$donor_id))

annot <- ensembldb::select(org.Hs.eg.db, keys = rownames(count), keytype = "ENSEMBL", columns = c("SYMBOL","ENSEMBL","ENTREZID"))

print(dim(annot[is.na(annot$SYMBOL),]))
print(dim(annot[duplicated(annot$SYMBOL),]))

annot <- annot[match(rownames(count),annot$ENSEMBL),]
print(head(annot))

mito_rename_dict <- c(
"ND1" = "MT-ND1", "ND2" = "MT-ND2", "COX1" = "MT-CO1", "COX2" = "MT-CO2",
"ATP8" = "MT-ATP8", "ATP6" = "MT-ATP6", "COX3" = "MT-CO3", "ND3" = "MT-ND3",
"ND4L" = "MT-ND4L", "ND4" = "MT-ND4", "ND5" = "MT-ND5", "ND6" = "MT-ND6",
"CYTB" = "MT-CYB", "RNR1" = "MT-RNR1", "RNR2" = "MT-RNR2",
"MT-ND1" = "MT-ND1", "MT-CO1" = "MT-CO1"
)

annot <- annot %>% mutate(rowname = case_when(
	SYMBOL %in% names(mito_rename_dict) ~ mito_rename_dict[SYMBOL],
	is.na(SYMBOL) ~ ENSEMBL,
	TRUE ~ SYMBOL) %>% unname())
rownames(count) <- make.names(annot$rowname, unique = T)
print(count[1:5,1:5])
export(annot, snakemake@output[["annot"]])

## data from paper (https://www-pnas-org.libproxy1.nus.edu.sg/doi/10.1073/pnas.1907883116), Figure S1
meta2 <- data.frame(donor_id = c(paste0("SC", 1:7), paste0("CT", 1:5)), sample_id =  c(paste0("SC", 1:7), paste0("CT", 1:5)),
			age = c(rep("110s", 7), "50s", "70s", "60s", "70s", "80s"),
			sex = c("female","male", rep("female", 4), "male", "male", rep("female", 3), "male"),
			self_reported_ethnicity = "unknown")
print(meta2)


metadata <- merge(metadata, meta2, by = "donor_id", all.x = T)
metadata <- metadata %>% column_to_rownames("cell_ID")
print(head(metadata))

seu <- CreateSeuratObject(counts = count, meta.data = metadata)
seu$percent_mito <- PercentageFeatureSet(object = seu, pattern = "^MT-")


print(head(seu@meta.data))
print(rownames(seu))
print(table(seu$donor_id, useNA = "always"))
print(dim(metadata[!is.na(metadata$donor_id),]))

seu <- subset(x = seu, subset = !is.na(donor_id))
print(table(seu$donor_id, useNA = "always"))

saveRDS(seu, snakemake@output[["seu"]])
export(seu@meta.data, snakemake@output[["meta"]])
