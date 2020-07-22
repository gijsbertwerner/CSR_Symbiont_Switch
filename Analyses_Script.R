#Gijsbert Werner, University of Oxford
#July 22, 2020


# Loading packages --------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ape)

# Loading data and trees --------------------------------------------------

dat_CSR_symb <-
  read.csv(file = "./Data/Plant species with symbiotic type and CSR strategy_Cosme et al 09-07-2020.csv",
           as.is = T,
           strip.white = T)
head(dat_CSR_symb)
sapply(dat_CSR_symb, class)

zanne_tree<-read.tree("./Data/Vascular_Plants_rooted.dated.tre")
zanne_tree

smith_brown_tree<-read.tree("./Data/ALLMB.tre")
smith_brown_tree
