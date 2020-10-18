# script to read and parse data files
rm(list=ls())
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

#############################################################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#############################################################
chromosome <- "chr21"
binsize <- 10000  #100000  #10000

#############################################################
hap_full_file <- "./../example_data/hap_full_scaffold_nov10_RPE1_22_"
het_file <- "./../example_data/het_coverage_jan3_RPE1_illumina_"

#############################################################
file_string <- paste(hap_full_file,chromosome,".dat",sep="")
hap_data2 <- read_delim(file_string ,delim="\t",col_names = FALSE)
raw_data2 <- dplyr::rename(hap_data2,"chr"="X1","arm"="X2","pos"="X3","ref_hash"="X4","alt_hash"="X5","hap"="X6","bl_string"="X8","block_hap"="X9","block"="X10")
hdata2 <- raw_data2 %>% separate(alt_hash,c("p","rbase","vbase"),"_",extra = "drop") %>% select("pos","rbase","vbase","hap","block")
hdata2 <- hdata2 %>% mutate(bin = floor(pos/binsize))

#############################################################
file_string <- paste(het_file,chromosome,".dat",sep="")
het_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
raw_data <- dplyr::rename(het_data,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
cdata <- raw_data %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
cdata <- cdata %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

jdata <- left_join(hdata2,cdata,by = c("pos","pos"))
jdata <- jdata %>% mutate(hapA_base = case_when(hap==1 ~ rbase, hap==-1 ~ vbase))
jdata <- jdata %>% mutate(hapB_base = case_when(hap==1 ~ vbase, hap==-1 ~ rbase))
jdata <- jdata %>% mutate(hapA_cov = case_when(hapA_base=="A" ~ as.numeric(abase), hapA_base=="C" ~ as.numeric(cbase), hapA_base=="G" ~ as.numeric(gbase),  hapA_base=="T" ~ as.numeric(tbase) ))
jdata <- jdata %>% mutate(hapB_cov = case_when(hapB_base=="A" ~ as.numeric(abase), hapB_base=="C" ~ as.numeric(cbase), hapB_base=="G" ~ as.numeric(gbase),  hapB_base=="T" ~ as.numeric(tbase) ))
jdata <- drop_na(jdata)

binned <- jdata %>% group_by(bin) %>% dplyr::summarise(A = mean(hapA_cov), B = mean(hapB_cov), U = mean(tot_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
binned <- binned %>% mutate(middle_pos = ((min_pos + max_pos)/2.0)) %>% filter(n>5)
binned_gather <- binned %>% gather(A,B,U,key = "hap",value = "cov")
binned_gather <- binned_gather %>% mutate(copy_number = cov/5.5)
binned_gather$hap = factor(binned_gather$hap, levels=c('U','A','B'))   # levels=c('A','B','U')

#############################################################
chromosome1_max = max(binned$max_pos)
ylimits <- c(0,3.5)   #c(0,50)  #c(0,180)
xlimits <- c(0,chromosome1_max)

#xlab <- c( 0, 30, 60, 90, 120, 150, 180, 210, 240)
xlab <- c( 0, 10, 20, 30, 40)

#############################################################
p1 <- ggplot(data=binned_gather) + geom_point(aes(x = middle_pos,y = copy_number,color = hap),size=0.5) +
  scale_color_manual(values=c("grey35", "firebrick2", "dodgerblue3")) +
  #facet_grid(.~hap) +
  facet_grid(hap~.) +
  scale_y_continuous(limits = ylimits,expand = c(0.0,0.0)) +
  scale_x_continuous(limits = xlimits,labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0)) +
  xlab(chromosome) +
  ylab("Copy Number") +
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"),legend.position = "none",axis.text.x = element_text(angle = 90),panel.grid.major = element_line(colour="grey", size=0.5))


p1

#############################################################
#chromosome1_scale = chromosome1_max/5000000
#output_file <- paste("~/2019_11_november_workdir/HCC1954/phased_hic_nov/phased_cn_",chromosome,".pdf",sep = "")
#ggsave(output_file, width = chromosome1_scale, height = 5, units = "cm")
