#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library(grid)
library(data.table)
library(cowplot)
library(patchwork)

# check is rstudio and if so set the working direcotry to curdir 
if(rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

df= fread("nobackups/chm13_clint.paf.tbl"); df$cutid = as.factor(df$cutid)


# make the identity histogram
bins = min(250, nrow(df)/10)
p1 = ggplot(data=df, aes(identity, fill=cutid, weight=aln_len))+ 
  geom_histogram(bins=bins) +
  ylab("bp") +
  scale_fill_brewer(palette = "Spectral", direction=-1) +
  scale_color_brewer(palette = "Spectral", direction=-1) +
  theme_classic()  + theme(legend.position = "none") +
  scale_y_continuous(label=comma)

# make the dot plot 
p2 = ggplot(data=df, aes(x=x1, xend=x2, y=y1, yend=y2, fill=cutid, color=cutid)) + 
  geom_segment() +   coord_fixed(ratio = 1) +
  scale_fill_brewer(palette = "Spectral", direction=-1) +
  scale_color_brewer(palette = "Spectral", direction=-1) +
  theme_classic()  + theme(legend.position = "none") +
  ylab(unique(df$q_name)) + xlab(unique(df$t_name)) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma)

# combine them 
p3 = p1/p2 + plot_layout(heights = c(1, 3)); p3

ggsave("~/Desktop/cen.png", plot=p3, height = 12, width = 8, dpi = 400)
