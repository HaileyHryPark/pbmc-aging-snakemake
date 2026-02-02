print(.libPaths())
.libPaths("resources/r_package")
print(.libPaths())

library(rio)
library(dplyr)
library(tidyverse)
library("DEswan")

data <- import(snakemake@input[["data"]])
f_data <- data %>% filter(sex == "Female")
m_data <- data %>% filter(sex == "Male")

print(f_data[1:15,1:15])
print(m_data[1:15,1:15])
print(unique(data$sex))

data_list <- list(f_data, m_data, data)
names(data_list) <- c("Female", "Male", "Both")

res <- lapply(as.list(names(data_list)), function(g){

	d <- data_list[[g]]
	
	start_time <- Sys.time()

	res.DEswan = DEswan(data.df = d[,-c(1:4)],
		qt = d$age,
		window.center = seq(20, 70, 1),
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
