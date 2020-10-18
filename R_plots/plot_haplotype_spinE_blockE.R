# script to read and parse data files
rm(list=ls())
library(readr)
library(tidyr)
library(dplyr)
library(matlab)
library(ggplot2)

###################################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

###################################
#chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
chrom_list <- c("chr21")
hap_file <- "./../example_data/hap_solution_nov7_RPE1_"

###################################
load_tenx_solution <- function(chrom_list,hap_file) {
  hap_tracks = tibble()
  for (chr in chrom_list) {
    file_string <- paste( hap_file , chr , ".dat" , sep="" )
    hap_data <- read_delim( file_string , delim="\t", col_names = FALSE )
    raw_data <- dplyr::rename(hap_data,"index"="X1","pos"="X2","ref_het_str"="X3","var_het_str"="X4","hap"="X5","rcount"="X6","vcount"="X7","spinE"="X8","blockE"="X9","block"="X10")
    raw_data <- raw_data %>% mutate(chromosome = chr)
    hap_tracks <- bind_rows(hap_tracks,raw_data)
  }
  return(hap_tracks)
}

########################################
########################################
break_list_linspace = linspace(0, 25000, n = 2000)

########################################
hap_tracks <- load_tenx_solution(chrom_list,hap_file)

########################################
hap_tracks_blockE <- hap_tracks %>% transform(bin = cut(-blockE,breaks=break_list_linspace)) %>% group_by(bin) %>% summarise(count = n(),midpoint = (max(-blockE) + min(-blockE))/2.0) %>% drop_na()
hap_tracks_spinE <- hap_tracks %>% transform(bin = cut(-spinE,breaks=break_list_linspace)) %>% group_by(bin) %>% summarise(count = n(),midpoint = (max(-spinE) + min(-spinE))/2.0) %>% drop_na()

p1_spinE <- ggplot(data=hap_tracks_spinE) + geom_line(aes(x=midpoint,y=count)) + #+ geom_freqpoly(aes(x=-spinE,y = ..count..),size=1.0,bins=90) +  #,center = 0 #pad = FALSE  #,na_rm=TRUE,inherit.aes = FALSE
  #+ geom_histogram(aes(x=-spinE),bins=90) +
  #+ geom_freqpoly(aes(x=-spinE),size=1.1,bins=90) +
  scale_x_continuous(expand = c(0.03,0.03),limits=c(0,2000)) + #2000  #3000  #limits=c(0,3000)  #expand = c(0.03,0.03) #expand = c(0.0,0.0)
  scale_y_continuous(expand = c(0,0)) +
  ylab("Count") +
  xlab("spinE" ) +
  theme_classic() +
  theme(legend.position="none")

p2_blockE <- ggplot(data=hap_tracks_blockE) + geom_line(aes(x=midpoint,y=count)) + #+ geom_freqpoly(aes(x=-blockE,y = ..count..),size=1.0,bins=90) + # ,binwidth = 1/4 #+ geom_histogram(aes(x=-blockE),bins=90) +  #,na_rm=TRUE,inherit.aes = FALSE
  #+ geom_histogram(aes(x=-blockE),bins=90) +
  #+ geom_freqpoly(aes(x=-blockE),size=1.1,bins=90) +
  geom_vline(aes(xintercept=700),size=0.8,alpha = 0.8,linetype = "dotted") + #,color="blue" #"700" #geom_vline(aes(xintercept=100,color="red"),size=1.1,alpha = 0.5) + #"100"   #geom_vline(aes(xintercept=2000,color="green"),size=1.1,alpha = 0.5) + #"2000"
  scale_x_continuous(expand = c(0.03,0.03),limits=c(0,9000)) +
  ylab("Count") +
  xlab("blockE") +
  theme_classic() +
  theme(legend.position="none")

#p1_spinE
p2_blockE
