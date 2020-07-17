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

H = 12
W = 8

df= fread("nobackups/gl8v8.tbl.tbl"); df$cutid = as.factor(df$cutid)


# make the identity histogram
bins = min(250, nrow(df)/10)
p1 = ggplot(data=df, aes(identity, fill=cutid, weight=aln_len))+ 
  geom_histogram(bins=bins) +
  ylab("bp") +
  scale_fill_brewer(palette = "Spectral", direction=-1) +
  scale_color_brewer(palette = "Spectral", direction=-1) +
  theme_cowplot()  + theme(legend.position = "none") +
  scale_y_continuous(label=comma)


small = df
nalns = 10e6
if(nrow(small) > nalns){
  small = small[sample(nrow(small), nalns),] # makes things too slow #, prob=small$matches/sum(small$matches)), ]
}
size = 0.5*(25.4 * W)/(nrow(small)); size


# make the dot plot 
p2 = ggplot(data=small, aes(x=x1, xend=x2, y=y1, yend=y2, fill=cutid, color=cutid)) + 
  geom_segment() +   coord_fixed(ratio = 1) +
  scale_fill_brewer(palette = "Spectral", direction=-1) +
  scale_color_brewer(palette = "Spectral", direction=-1) +
  theme_cowplot()  + theme(legend.position = "none") +
  ylab(unique(df$q_name)) + xlab(unique(df$t_name)) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma)

# combine them 
p3 = p1/p2 + plot_layout(heights = c(1, 3)); p3


density= 50
ggsave("~/Desktop/example.png", plot=p3, height = density*H, width = density*W, units="mm", limitsize = F)


