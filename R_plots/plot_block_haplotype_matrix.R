# script to read and parse data files
rm(list=ls())
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

########################################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

########################################
bmatrix_file <- "./../example_data/block_phase_matrix_nov10_RPE1_22_chr21.dat"

########################################
bmatrix_data <- read_delim( bmatrix_file , delim="\t", col_names = FALSE )
bmatrix_data <- dplyr::rename( bmatrix_data ,"i"="X1","j"="X2","iraw"="X3","jraw"="X4","hap_i"="X5","hap_j"="X6","Bij"="X7")
bmatrix_data <- bmatrix_data %>% mutate(phasedij=Bij*hap_i*hap_j)
bmatrix_data <- bmatrix_data %>% mutate(sign_phasedij = sign(phasedij))

p1 <- ggplot(bmatrix_data, aes(x = i, y = j)) +
  geom_tile(aes(fill=sign_phasedij)) +
  scale_fill_gradientn(colours=c("blue","gray80","red"), limits=c(-1,1)) +
  labs(fill = "sign(Mij)") +
  theme_classic() 

p1
