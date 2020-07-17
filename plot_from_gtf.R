#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
#library(dplyr)
library(grid)
#library(gridBase)
library(gridExtra)
library(data.table)
library(gtable)
library(tidyr)
library(dplyr)



names = c('contig', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
df = fread('~/Desktop/EichlerVolumes/chm13_t2t/nobackups/tmp.gff', col.names = names);df

df = data.table(df %>% 
  extract(attributes, "gene_id", 'transcript_id=([^=;]+);',remove = FALSE) %>% 
  extract(attributes, "id", 'gene_id=([^=;]+);',remove = FALSE) %>% 
  extract(attributes, "name", 'gene_name=([^=;]+);',remove = FALSE))
df$gene_id=as.factor(df$gene_id)

d = df %>% filter(feature == "transcript") %>%
  mutate(range = end-start) %>%
  group_by(id) %>%
  filter(range == max(range)) %>%
  ungroup() %>%
  distinct(gene_id)
  # %>% summarise(x=max( max(end)-min(start)))
df = df %>% filter(gene_id %in% d$gene_id)

p=ggplot()+
  geom_segment(data = df[df$feature=='exon',], aes(x=start, xend=end, y=gene_id, yend=gene_id), size=4)+
  geom_segment(data = df[df$feature=='CDS',], aes(x=start, xend=end, y=gene_id, yend=gene_id), size=8)+
  geom_segment(data = df[df$feature=='transcript',], aes(x=start, xend=end, y=gene_id, yend=gene_id), size=1)+
  scale_x_continuous( expand = c(0, 0), limits = c(1,max(df$end)),breaks = c(1,max(df$end)) ) +
  ylab('')+xlab('')+
  theme_classic()+  
  theme(plot.margin=margin(1,1,1,1, 'cm'), axis.line.y =  element_blank(), axis.ticks.y = element_blank());p


h=min(6,length(unique(df$gene_id)))
ggsave('~/Desktop/gene.pdf', plot = p, width = 8, height = h, units = 'in')
