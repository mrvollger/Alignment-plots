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

H = 16
W = 12

df= fread("all.cen.paf.tbl"); df$cutid = as.factor(df$cutid)
#df$cutid = factor(df$cutid, levels = seq(9,0,-1))
direction = -1
df=df[order(df$cutid)]

names = unique(c(df$t_name, df$q_name))
df$t_name=factor(df$t_name, levels=names)
df$q_name=factor(df$q_name, levels=names)

if(F){
temp = sample_n(df,100)
temp$x1=0
temp$x2=5000
temp$y1=0
temp$y2=5000
temp=temp[order(temp$cutid),]


ggplot(data=temp, aes(x=x1, xend=x2, y=y1, yend=y2, fill=cutid, color=cutid)) + 
  geom_segment() +   coord_fixed(ratio = 1) +
  scale_fill_brewer(palette = "Spectral", direction=direction) +
  scale_color_brewer(palette = "Spectral", direction=direction) +
  theme_cowplot()  + theme(legend.position = "none") +
  ylab(unique(df$q_name)) + xlab(unique(df$t_name)) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma)
}


# make the identity histogram
bins = min(50, max(nrow(df)/10,10))
p1 = ggplot(data=df, aes(identity, fill=cutid, weight=aln_len))+ 
  geom_histogram(bins=bins) +
  ylab("bp") +
  scale_fill_brewer(palette = "Spectral", direction=direction) +
  scale_color_brewer(palette = "Spectral", direction=direction) +
  theme_cowplot()  + theme(legend.position = "none") +
  scale_y_continuous(label=comma)
p1

small = df
nalns = 10e6
if(nrow(small) > nalns){
  small = small[sample(nrow(small), nalns),] # makes things too slow #, prob=small$matches/sum(small$matches)), ]
}
size = 0.5*(25.4 * W)/(nrow(small)); size
small=small[order(small$cutid),]


# make the dot plot 
p2 = ggplot(data=small, aes(x=x1, xend=x2, y=y1, yend=y2, fill=cutid, color=cutid)) + 
  geom_segment(size=1) +   coord_fixed(ratio = 1) +
  scale_fill_brewer(palette = "Spectral", direction=direction) +
  scale_color_brewer(palette = "Spectral", direction=direction) +
  theme_bw()  +
  theme(legend.position = "none") +
  #ylab(unique(df$q_name)) + xlab(unique(df$t_name)) +
  ylab("") + xlab("") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma)+facet_grid(vars(q_name), vars(t_name))

# combine them 
p3 = p1/p2 + plot_layout(heights = c(1, 4));

p3

density= 1
ggsave("~/Desktop/diag_all.pdf", plot=p3, height = density*H, width = density*W, limitsize = F)







