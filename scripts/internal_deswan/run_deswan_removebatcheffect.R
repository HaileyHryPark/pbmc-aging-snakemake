print(.libPaths())
.libPaths("resources/r_package")
print(.libPaths())


library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(limma)
library("DEswan")

data <- import(snakemake@input[["data"]])
f_data <- data %>% filter(sex == "female")
m_data <- data %>% filter(sex == "male")


limma_data_list <- lapply(list(f_data, m_data, data), function(d){

limma_data <- d %>% select(-c(rowname,age,sex,dataset,ethnicity))
limma_res <- removeBatchEffect(t(limma_data), d$dataset)
limma_res <- as.data.frame(t(limma_res))
limma_res$rowname <- d$rowname
limma_res$age <- d$age
limma_res$sex <- d$sex
limma_res$dataset <- d$dataset
limma_res$ethnicity <- d$ethnicity

return(limma_res)

})

names(limma_data_list) <- c("Female", "Male", "Both")

res <- lapply(as.list(names(limma_data_list)), function(g){

	d <- limma_data_list[[g]]

        start_time <- Sys.time()

        res.DEswan= DEswan(data.df = d %>% select(-c(rowname,age,sex,dataset,ethnicity)),
                qt = d$age,
                window.center = seq(20, 90, 1),
                buckets.size = 20)

        end_time <- Sys.time()
        print(end_time-start_time)
        print(head(res.DEswan$p))
        print(head(res.DEswan$coeff))

	res.DEswan.wide.coeff=reshape.DEswan(res.DEswan,parameter = 2,factor = "qt") %>% as.data.frame() %>% mutate(gender = g)
	res.DEswan.wide.p=reshape.DEswan(res.DEswan,parameter = 1,factor = "qt") %>% as.data.frame() %>% mutate(gender = g)
        res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH") %>% as.data.frame() %>% mutate(gender = g)

        return(list(res.DEswan, res.DEswan.wide.coeff, res.DEswan.wide.p, res.DEswan.wide.q))

})

all_res <- lapply(res, `[[`, 1)
names(all_res) <- c("Female", "Male", "Both")
saveRDS(all_res, snakemake@output[["res"]])

export(bind_rows(lapply(res, `[[`, 2)), snakemake@output[["coef"]])
export(bind_rows(lapply(res, `[[`, 3)), snakemake@output[["p"]])
export(bind_rows(lapply(res, `[[`, 4)), snakemake@output[["q"]])
