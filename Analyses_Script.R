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
nrow(dat_CSR_symb)

zanne_tree <- read.tree("./Data/Vascular_Plants_rooted.dated.tre")
zanne_tree

smith_brown_tree <- read.tree("./Data/ALLMB.tre")
smith_brown_tree

# Data Cleaning -----------------------------------------------------------

#Species formatting
head(dat_CSR_symb$Species_name)
head(zanne_tree$tip.label)
head(smith_brown_tree$tip.label)
dat_CSR_symb$Species_name <- gsub(pattern = " ",
                                  replacement = "_",
                                  dat_CSR_symb$Species_name)

#How many of the 3014 species in the database are absent in Zanne?
length(setdiff(dat_CSR_symb$Species_name, zanne_tree$tip.label))

#How many of the 3014 species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb$Species_name, smith_brown_tree$tip.label))
#Only 320 of 3014 lacking. I.e. ~90% is present. That seems good enough for now.
write.csv(
  setdiff(dat_CSR_symb$Species_name, smith_brown_tree$tip.label),
  quote = F,
  row.names = F,
  file = "./Data/Cosme_Species_Missing_Smith_Tree.csv"
)

#Extract the appropriate subtree from Smith&Brown
analysis_tree <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label,dat_CSR_symb$Species_name)
  )
analysis_tree
