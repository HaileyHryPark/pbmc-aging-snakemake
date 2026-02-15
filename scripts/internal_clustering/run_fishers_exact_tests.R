library(rio)
library(tidyverse)


clust_f <- import(snakemake@input[["clust_f"]])
clust_m <- import(snakemake@input[["clust_m"]])
gs <- readRDS(snakemake@input[["gs"]])
ias_gs <- gs$Serrano_Iron_Accumulation_Geneset
smy_gs <- gs$SenMaYO

sink(snakemake@output[["res1"]])
print("### IAS")
# Gene level
print("Gene level")
a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique(), ias_gs))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique()) - a
c <- length(ias_gs)- a
d <- length(clust_f %>% pull(gene) %>% unique()) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

# Feature level
print("Feature level")
ias_features <- clust_f %>% filter(gene %in% ias_gs) %>% pull(feature)
print(ias_features)

a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature), ias_features))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature)) - a
c <- length(ias_features) - a
d <- nrow(clust_f) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

print("### SenMaYO")
# Gene level
print("Gene level")
a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique(), smy_gs))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique()) - a
c <- length(smy_gs)- a
d <- length(clust_f %>% pull(gene) %>% unique()) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

# Feature level
print("Feature level")
smy_features <- clust_f %>% filter(gene %in% smy_gs) %>% pull(feature)
print(smy_features)

a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature), smy_features))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature)) - a
c <- length(smy_features) - a
d <- nrow(clust_f) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

sink()

sink(snakemake@output[["res2"]])
print("### IAS")
# Gene level
print("Gene level")
a <- length(intersect(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(gene) %>% unique(), ias_gs))
b <-  length(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(gene) %>% unique()) - a
c <- length(ias_gs)- a
d <- length(clust_f %>% pull(gene) %>% unique()) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

# Feature level
print("Feature level")
ias_features <- clust_f %>% filter(gene %in% ias_gs) %>% pull(feature)
print(ias_features)

a <- length(intersect(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(feature), ias_features))
b <-  length(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(feature)) - a
c <- length(ias_features) - a
d <- nrow(clust_f) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

print("### SenMaYO")
# Gene level
print("Gene level")
a <- length(intersect(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(gene) %>% unique(), smy_gs))
b <-  length(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(gene) %>% unique()) - a
c <- length(smy_gs)- a
d <- length(clust_f %>% pull(gene) %>% unique()) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

# Feature level
print("Feature level")
smy_features <- clust_f %>% filter(gene %in% smy_gs) %>% pull(feature)
print(smy_features)

a <- length(intersect(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(feature), smy_features))
b <-  length(clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(feature)) - a
c <- length(smy_features) - a
d <- nrow(clust_f) - a - b - c

mat <- matrix(c(a, c, b, d), nrow = 2)
print(mat)
fisher.test(mat, alternative = "greater")

sink()

## Continuous increase in female
sink(snakemake@output[["res3"]])
print("### Female continuous increase cluster")
print("#FCI_MEI")
all_features_n <- length(unique(c(clust_f$feature, clust_m$feature)))
fci <- clust_f %>% filter(final_cluster == "Continuous\nincrease") %>% pull(feature)
mei <- clust_m %>% filter(final_cluster == "Early\nincrease") %>% pull(feature)
fcimei <- intersect(fci, mei)

a <- length(fcimei)
b <- length(fci) - a
c <- length(mei) - a
d <- all_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Female = c("In_ci_cluster", "Not_in_ci_cluster"),
                Male   = c("Early_increase", "Not_early_increase")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("#CD8 FCI")

all_f_features_n <- nrow(clust_f)
fci
cd8_f <- clust_f %>% filter(celltype == "CD8 T") %>% pull(feature)
cd8_fci <- clust_f %>% filter(celltype == "CD8 T", final_cluster == "Continuous\nincrease") %>% pull(feature)

a <- length(cd8_fci)
b <- length(fci) - a
c <- length(cd8_f) - a
d <- all_f_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("CD8T", "Not_CD8T")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("#CD8 FCIMEI")

all_features_n
fcimei
cd8 <- unique(clust_f %>% filter(celltype == "CD8 T") %>% pull(feature), 
		clust_m %>% filter(celltype == "CD8 T") %>% pull(feature))
cd8_fcimei <- intersect(cd8, fcimei)

a <- length(cd8_fcimei)
b <- length(fcimei) - a
c <- length(cd8) - a
d <- all_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("CD8T", "Not_CD8T")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("### Female late increase cluster")

print("#T FLI")

all_f_features_n 
fli <- clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature)
t_f <- clust_f %>% filter(celltype %in% c("CD8 T", "CD4 T")) %>% pull(feature)
t_fli <- intersect(fli, t_f)

a <- length(t_fli)
b <- length(fli) - a
c <- length(t_f) - a
d <- all_f_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("T", "Not_T")
              ))

print(mat)

fisher.test(mat, alternative = "less")

print("#No male SAF FLI")

all_features_n
fli
only_f <- setdiff(clust_f$feature, clust_m$feature) 
only_f_fli <- intersect(fli, only_f)

a <- length(only_f_fli)
b <- length(fli) - a
c <- length(only_f) - a
d <- all_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("only_f", "Not_only_f")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("###Female inverted U-shape cluster")

print("#CD4 T FIU")

all_f_features_n 
fiu <- clust_f %>% filter(final_cluster == "Inverted\nUshape") %>% pull(feature)
cd4t_f <- clust_f %>% filter(celltype == "CD4 T") %>% pull(feature)
cd4t_fiu <- intersect(fiu, cd4t_f)

a <- length(cd4t_fiu)
b <- length(fiu) - a
c <- length(cd4t_f) - a
d <- all_f_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("CD4T", "Not_CD4T")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("#No male SAF FIU")

all_features_n
fiu
only_f
only_f_fiu <- intersect(fiu, only_f)

a <- length(only_f_fiu)
b <- length(fiu) - a
c <- length(only_f) - a
d <- all_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("only_f", "Not_only_f")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("###Male early fluctuation cluster")

print("#CD4 T MEF")

all_m_features_n <- nrow(clust_m)
mef <- clust_m %>% filter(final_cluster == "Early\nfluctuation") %>% pull(feature)
cd4t_m <- clust_m %>% filter(celltype == "CD4 T") %>% pull(feature)
cd4t_mef <- intersect(mef, cd4t_m)

a <- length(cd4t_mef)
b <- length(mef) - a
c <- length(cd4t_m) - a
d <- all_m_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("CD4T", "Not_CD4T")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

print("#No female SAF MEF")

all_features_n
mef
only_m <- setdiff(clust_m$feature, clust_f$feature)
only_m_mef <- intersect(mef, only_m)

a <- length(only_m_mef)
b <- length(mef) - a
c <- length(only_m) - a
d <- all_features_n - (a + b + c)

mat <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(
                Cluster = c("In_cluster", "Not_in_cluster"),
                CellType = c("only_m", "Not_only_m")
              ))

print(mat)

fisher.test(mat, alternative = "greater")

sink()



