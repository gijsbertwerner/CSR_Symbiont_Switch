#Gijsbert Werner, University of Oxford
#July 22, 2020

# Loading packages --------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ape)
library(diversitree)
library(qpcR)
library(corHMM)
library(phytools)
library(phylolm)
library(parallel)

# Loading data and trees --------------------------------------------------

dat_CSR_symb <-
  read.csv(file = "./Data/Plant species with symbiotic type and CSR strategy_Cosme et al 09-07-2020.csv",
           as.is = T,
           strip.white = T)
head(dat_CSR_symb)
dat_CSR_symb$Symbiotic_type <-
  gsub(pattern = "NM-AM",
       replacement = "NMAM",
       dat_CSR_symb$Symbiotic_type)
dat_CSR_symb$Symbiotic_type <-
  gsub(pattern = "EcM-AM",
       replacement = "EcMAM",
       dat_CSR_symb$Symbiotic_type)
dat_CSR_symb$Symbiotic_type <-
  gsub(pattern = "AM-Nod",
       replacement = "AMNod",
       dat_CSR_symb$Symbiotic_type)
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
#What we could still do to get the numbers up
# (1) Fuzzy matching of species names to tree -> so small spelling variations means it doesn't immediately drop out.
# (2) Manually check the missing 10% with reference to the tree. See if synonyms are present.
# @Marco: do you think either of these is worth it?

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
analysis_dat_CSR_symb[duplicated(analysis_dat_CSR_symb$Species_name),]
#Yes, some duplicates, but only about ~30. Small variations in CSR-values, agreement on symbiont type

#For now, simply remove the duplicates. Other option, take average values for CSR values.
analysis_dat_CSR_symb <-
  analysis_dat_CSR_symb[!duplicated(analysis_dat_CSR_symb$Species_name), ]
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

# Basic statistics on the C, S and R values
summary(analysis_dat_CSR_symb$C.selection)
summary(analysis_dat_CSR_symb$S.selection)
summary(analysis_dat_CSR_symb$R.selection)

ggplot(data = analysis_dat_CSR_symb) +
  geom_histogram(aes(C.selection))

ggplot(data = analysis_dat_CSR_symb) +
  geom_freqpoly(aes(C.selection), colour = "Red") +
  geom_freqpoly(aes(S.selection), colour = "Blue") +
  geom_freqpoly(aes(R.selection), colour = "Green")
#Ok, so looking at the graphs quite a lot of plants on the 0 side for particularly S and R selection.

#Ways to treat CSR analytically:
#1. Turn it into a categorical variable: assign one of three categorical states based on what it is most selected for.
#2. Treat as three distinct continuous variables, and repeat the analyses for each.
#Neither of these captures perfectly what we are trying to measure, but they may do for our purposes. Let's have a look.

# Analyses ----------------------------------------------------------------

#Let's start with ASR (Ancestral State Reconstructions) for symbiont state, and plot them on the tree

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type)

#Run ASRs
ASR_symnbiont_type_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 100,
    n.cores = 6
  )
ASR_symnbiont_type_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 100,
    n.cores = 6
  )
ASR_symnbiont_type_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 100,
    n.cores = 6
  )

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symnbiont_type_ER_yang$AICc,
    ASR_symnbiont_type_ARD_yang$AICc,
    ASR_symnbiont_type_SYM_yang$AICc
  )
)

#ARD by far the best
ASR_symnbiont_type_ARD_yang
#Consideration: we have lots of states (8), and some have only few cases.
#We could simplify the inference by dropping some states and/or subsuming them in other states.
#That would make phylogenetic inference easier, and the results more interpretable, but also lose biological detail. 



#Steps
#Plot categories of symbiotic state on outline
#ASR of categorical states
#ASRs of the three CSR values. 3x
#Plot these on top of each other
#Pagel's model, need two binary variables.
#Turn csr into three categorical ones, see if correlated three variables with corhmm.
