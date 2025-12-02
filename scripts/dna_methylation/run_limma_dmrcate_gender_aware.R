library(rio)
library(DMRcate)
library(limma)
library(tidyverse)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

### Functions
# Design matrix
make_design <- function(pheno, mode = c("both", "sexstrat")) {
  mode <- match.arg(mode)

  if (mode == "both") {
    # Formula: age_group * sex
    design <- model.matrix(~ AgeGroup * Gender, data = pheno)
  } else {
    # Single gender: ~ age_group
    design <- model.matrix(~ AgeGroup, data = pheno)
  }
  return(design)
}

# Define contrasts
contrast_age <- function(design, mode = c("both", "sexstrat")) {
  mode <- match.arg(mode) 

  print(colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))
  print(colnames(design))

  c1 <- makeContrasts(AG_40to60_vs_lt40 = `AgeGroup40to60`,
                      levels = design)
  print(c1)

  c2 <- makeContrasts(AG_gt60_vs_40to60 = `AgeGroupgt60` - `AgeGroup40to60`,
                      levels = design)

  c3 <- makeContrasts(AG_gt60_vs_lt40 = `AgeGroupgt60`,
                      levels = design)

  if(mode == "both"){
    c4 <- makeContrasts(Interaction_40to60 = `AgeGroup40to60_Gendermale`,
                        levels = design)

    c5 <- makeContrasts(Interaction_gt60 = `AgeGroupgt60_Gendermale`,
                        levels = design)

    c6 <- makeContrasts(Interaction_40to60_gt60 = `AgeGroupgt60_Gendermale` - `AgeGroup40to60_Gendermale`,
                        levels = design)

    return(list(AG_40to60_vs_lt40 = c1, AG_gt60_vs_40to60 = c2, AG_gt60_vs_lt40 = c3, Interaction_40to60 = c4, Interaction_gt60 = c5, Interaction_40to60_gt60 = c6))
  }else{
    list(AG_40to60_vs_lt40 = c1, AG_gt60_vs_40to60 = c2, AG_gt60_vs_lt40 = c3)
  }
}

# Running limma & DMRcate
run_pipeline <- function(M, pheno, design_mode) {

  design <- make_design(pheno, mode = design_mode)
  ctr <- contrast_age(design, mode = design_mode)

  limma_list <- list()
  dmr_list <- list()

  for (contrast_name in names(ctr)) {
    cat("\n### Running contrast:", contrast_name, "###\n")

    fit <- lmFit(M, design)
    fit2 <- contrasts.fit(fit, ctr[[contrast_name]])
    fit2 <- eBayes(fit2)

    limma_res <- topTable(fit2, coef = 1, number = Inf) %>% 
	rownames_to_column("probeID") %>% 
	left_join(annot, by = join_by("probeID" == "Name"))
    limma_list[[contrast_name]] <- limma_res

    ann <- cpg.annotate(
      datatype = "array",
      object = M,
      what = "M",
      arraytype = "450K",
      analysis.type = "differential",
      design = design,
      contrasts = TRUE,
      cont.matrix = ctr[[contrast_name]],
      coef = contrast_name
    )

    dmrs <- dmrcate(ann, lambda = 1000, C = 2)

    dmr_list[[contrast_name]] <- extractRanges(dmrs)
  }

  return(list(limma = limma_list, dmr = dmr_list))
}

### Main
# Read data
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame()
print(head(annot))

data <- import(snakemake@input[["mvalue"]])
metadata <- import(snakemake@input[["metadata"]]) %>% filter(!is.na(Age), Gender != "")
metadata$AgeGroup <- factor(metadata$AgeGroup, levels = c("<40","40-60",">60"), labels = c("lt40","40to60","gt60"))
metadata$Gender <- factor(metadata$Gender, levels = c("female", "male"))

mat <- data %>% column_to_rownames("probe_prefix") %>% as.matrix()

print(head(mat))
print(dim(mat))
print(dim(metadata))

print(setdiff(rownames(mat), annot$Name))

# Subset metadata to match NPX matrix
meta_b <- metadata
mat_b <- mat[, meta_b$Basename]

meta_f <- metadata %>% filter(Gender == "female")
mat_f <- mat[, meta_f$Basename]

meta_m <- metadata %>% filter(Gender == "male")
mat_m <- mat[, meta_m$Basename]

res_both <- run_pipeline(mat_b, meta_b, "both")

res_female <- run_pipeline(mat_f, meta_f, "sexstrat")

res_male <- run_pipeline(mat_m, meta_m, "sexstrat")

saveRDS(res_both,   snakemake@output[["res_both"]])
saveRDS(res_female, snakemake@output[["res_female"]])
saveRDS(res_male,   snakemake@output[["res_male"]])
