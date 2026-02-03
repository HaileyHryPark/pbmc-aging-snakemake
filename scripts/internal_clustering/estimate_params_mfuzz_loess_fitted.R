library(rio)
library(Mfuzz)
library(dplyr)
library(tidyverse)

eset <- table2eset(filename = snakemake@input[["mfuzz_mat"]])

set.seed(123)
m1 <- mestimate(eset)
print("M1:")
print(m1)

set.seed(123)
 plot <-
   Dmin(
     eset,
     m = m1,
     crange = seq(2, 40, 2),
     repeats = 3,
     visu = TRUE
   )

df <- plot %>% data.frame(distance = plot, k = seq(2, 40, 2))

 plot <- df %>% 
   ggplot(aes(k, distance)) +
   geom_point(shape = 21, size = 4, fill = "black") +
   # geom_smooth() +
   geom_segment(aes(
     x = k,
     y = 0,
     xend = k,
     yend = distance
   )) +
   theme_bw() +
   theme(
     # legend.position = c(0, 1),
     # legend.justification = c(0, 1),
     panel.grid = element_blank(),
     axis.title = element_text(size = 13),
     axis.text = element_text(size = 12),
     panel.background = element_rect(fill = "transparent", color = NA),
     plot.background = element_rect(fill = "transparent", color = NA),
     legend.background = element_rect(fill = "transparent", color = NA)
   ) +
   labs(x = "Cluster number",
        y = "Min. centroid distance", 
	caption = paste("m1:", m1)) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

 plot

 ggsave(plot,
        filename = snakemake@output[["plot"]],
        width = 4,
        height = 4)
export(df, snakemake@output[["table"]])
