print(.libPaths())
.libPaths(c("/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package", "/home/users/nus/e0859928/opt/miniforge3/envs/rbase/lib/R/library", .libPaths()))
print(.libPaths())

library(rio)
library(dplyr)
library(tidyverse)
library("DEswan")

data <- import(snakemake@input[["data"]])

ds_res <- lapply(as.list(unique(data$dataset)), function(ds){

data <- data %>% filter(dataset == ds)
f_data <- data %>% filter(dataset == ds, sex == "female")
m_data <- data %>% filter(dataset == ds, sex == "male")

print(f_data[1:15,1:15])
print(m_data[1:15,1:15])
print(unique(data$sex))

data_list <- list(f_data, m_data, data)
names(data_list) <- c("Female", "Male", "Both")

res <- lapply(as.list(names(data_list)), function(g){

	d <- data_list[[g]]
	if(nrow(d) < 10) return(list(NULL, NULL, NULL, NULL))
	
	start_time <- Sys.time()

	res.DEswan = DEswan(data.df = d[,-c(1:5)],
		qt = d$age,
		window.center = seq(min(d$age), max(d$age), 1),
		buckets.size = 20)

	end_time <- Sys.time()
	print(end_time-start_time)
	print(head(res.DEswan$p))
	print(head(res.DEswan$coeff))

	res.DEswan.wide.coeff=reshape.DEswan(res.DEswan,parameter = 2,factor = "qt") %>% as.data.frame() %>% mutate(gender = g, dataset = ds)
	res.DEswan.wide.p=reshape.DEswan(res.DEswan,parameter = 1,factor = "qt") %>% as.data.frame() %>% mutate(gender = g, dataset = ds)
	res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH") %>% as.data.frame() %>% mutate(gender = g, dataset = ds)

	return(list(res.DEswan, res.DEswan.wide.coeff, res.DEswan.wide.p, res.DEswan.wide.q))
})

coef <- lapply(res, `[[`, 2)
coef <- coef[!unlist(lapply(coef, is.null))]

p <- lapply(res, `[[`, 3)
p <- p[!unlist(lapply(p, is.null))]

q <- lapply(res, `[[`, 4)
q <- q[!unlist(lapply(q, is.null))]

return(list(bind_rows(coef),bind_rows(p),bind_rows(q)))
})

export(bind_rows(lapply(ds_res, `[[`, 1)), snakemake@output[["coef"]])
export(bind_rows(lapply(ds_res, `[[`, 2)), snakemake@output[["p"]])
export(bind_rows(lapply(ds_res, `[[`, 3)), snakemake@output[["q"]])
