print(.libPaths())
.libPaths(c("/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package", "/home/users/nus/e0859928/opt/miniforge3/envs/rbase/lib/R/library", .libPaths()))
print(.libPaths())

library(rio)
library(dplyr)
library(tidyverse)
library("DEswan")

data <- import(snakemake@input[["data"]])
f_data <- data %>% filter(sex == "female")
m_data <- data %>% filter(sex == "male")

print(f_data[1:15,1:15])
print(m_data[1:15,1:15])
print(unique(data$sex))

data_list <- list(f_data, m_data, data)
names(data_list) <- c("Female", "Male", "Both")

res <- lapply(list(5,10,15,20,25), function(b){
	rres <- lapply(as.list(names(data_list)), function(g){

		d <- data_list[[g]]
	
		start_time <- Sys.time()

		res.DEswan = DEswan(data.df = d[,-c(1:5)],
			qt = d$age,
			window.center = seq(20, 90, 1),
			buckets.size = b, 
			covariates = d$dataset)

		end_time <- Sys.time()
		print(end_time-start_time)
		print(head(res.DEswan$p))
		print(head(res.DEswan$coeff))

		res.DEswan.wide.p=reshape.DEswan(res.DEswan,parameter = 1,factor = "qt")
		res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH") %>% as.data.frame() %>% mutate(gender = g, bucket = b)

		return(res.DEswan.wide.q)
	})
	return(bind_rows(rres))
})

all_res <- bind_rows(res) 
export(all_res, snakemake@output[["res"]])

