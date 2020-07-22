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
#About one third of the data set is not present in Zanne. That's a bit much.

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
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb$Species_name)
  )
analysis_tree

analysis_dat_CSR_symb <-
  dat_CSR_symb %>% filter(Species_name %in% analysis_tree$tip.label)
nrow(analysis_dat_CSR_symb)
#More species in data set, than in tree.
#This should not be possible. Could be due to duplicates in the dataset?
analysis_dat_CSR_symb[duplicated(analysis_dat_CSR_symb$Species_name), ]
#Yes, some duplicates, but only about ~30. Small variations in CSR-values, agreement on symbiont type

#For now, simply remove the duplicates. Other option, take average values for CSR values.
analysis_dat_CSR_symb <-
  analysis_dat_CSR_symb[!duplicated(analysis_dat_CSR_symb$Species_name),]
nrow(analysis_dat_CSR_symb)
#Now matches

#Lastly some sanity checks
table(analysis_dat_CSR_symb$Symbiotic_type) #This seems a reasonable distribution.
#Do they all sum to 100?
length(which(
  round(
    analysis_dat_CSR_symb$C.selection + analysis_dat_CSR_symb$S.selection +
      analysis_dat_CSR_symb$R.selection,
    0
  ) == 100
))


# Descriptives ------------------------------------------------------------

#Plot the 

# Basic statistics on the C, S and R values
summary(analysis_dat_CSR_symb$C.selection)
summary(analysis_dat_CSR_symb$S.selection)
summary(analysis_dat_CSR_symb$R.selection)

ggplot(data = analysis_dat_CSR_symb)+
  geom_histogram(aes(C.selection))

ggplot(data = analysis_dat_CSR_symb)+
  geom_freqpoly(aes(C.selection),colour="Red")+
  geom_freqpoly(aes(S.selection),colour="Blue")+
  geom_freqpoly(aes(R.selection),colour="Green")

#Steps
