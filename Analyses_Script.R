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
library(RColorBrewer)
library(Rphylopars)
library(viridis)

# Loading data and trees --------------------------------------------------

dat_CSR_symb <-
  read.csv(file = "./Data/Plant species with symbiotic type and CSR strategy_Cosme et al 16-09-2020.csv",
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
table(dat_CSR_symb$CSR_categorical_level)
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

#Print species from database not in tree.
write.csv(
  setdiff(dat_CSR_symb$Species_name, smith_brown_tree$tip.label),
  quote = F,
  row.names = F,
  file = "./Output/Cosme_Species_Missing_Smith_Tree.csv"
)
#Print all tree species
write.csv(
  smith_brown_tree$tip.label,
  quote = F,
  row.names = F,
  file = "./Output/SmithBrownAllSpecies.csv"
)

###Manual substitutions of data names
sub_table <- read.csv("./Data/SolvingMismatchesSmith.csv")
head(sub_table)

#Replace in the CSR data file
dat_CSR_symb$Species_name <-
  ifelse(
    is.na(
      match(dat_CSR_symb$Species_name, sub_table$Name.in.CSR.data.set)
    ),
    dat_CSR_symb$Species_name,
    sub_table$Equivalent.in.Smith.Brown.Tree[match(dat_CSR_symb$Species_name, sub_table$Name.in.CSR.data.set)]
  )
dat_CSR_symb$Species_name
length(setdiff(dat_CSR_symb$Species_name, smith_brown_tree$tip.label)) #Yes, now only 17 species missing. 


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

# Analyses ----------------------------------------------------------------

###Symbiont type

#Let's start with ASR (Ancestral State Reconstructions) for symbiont state, and plot them on the tree

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type)

#Run ASRs
ASR_symbiont_type_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_symbiont_type_ER_yang")
load("./Output/ASR_symbiont_type_ARD_yang")
load("./Output/ASR_symbiont_type_SYM_yang")

save(ASR_symbiont_type_ER_yang, file = "./Output/ASR_symbiont_type_ER_yang")
save(ASR_symbiont_type_ARD_yang, file = "./Output/ASR_symbiont_type_ARD_yang")
save(ASR_symbiont_type_SYM_yang, file = "./Output/ASR_symbiont_type_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_ER_yang$AICc,
    ASR_symbiont_type_ARD_yang$AICc,
    ASR_symbiont_type_SYM_yang$AICc
  )
)

#SYM by far the best
ASR_symbiont_type_SYM_yang
plotMKmodel(ASR_symbiont_type_SYM_yang)
table(analysis_dat_CSR_symb$Symbiotic_type) #States are numbered in the modeling: this is what types the numbers represent, they are ordered aphabetically, it sems.

#Consideration: we have lots of states (8), and some have only few cases.
#We could simplify the inference by dropping some states and/or subsuming them in other states.
#That would make inference easier, and the results more interpretable (not so many state transitions), but also lose biological detail.

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_symbiont_type <-
  analysis_dat_CSR_symb_ASR_symbiont_type %>%
  dplyr::select(Symbiotic_type)
row.names(dat_plot_symbiont_type) <-
  analysis_dat_CSR_symb_ASR_symbiont_type$Species_name
# dat_plot_symbiont_type$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_type$Symbiotic_type))
head(dat_plot_symbiont_type)

#Symbionts ASR - plot to pdf
pdf("./Output/ASRSymbiontType.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()

table(analysis_dat_CSR_symb$Symbiotic_type)
#Ok, what are we seeing here?
#A reconstruction of symbiont type: ancestral state is estimated as AM, and AM seems maintained througout.
#With transitions towards other types that are quite 'tippy' (i.e. recent/shallow) at this phylogenetic scale.
#No major suprises here, I would say. Everything in line with earlier work (including Nadia and mine previous paper. )

#######Explore ways of simplifying this ASR

#####First, lump AM+NMAM
#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped$Symbiotic_type<-
  ifelse(analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped$Symbiotic_type %in% c("AM","NMAM"),
         "AM_NMAM",analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped$Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped)
table(analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped$Symbiotic_type)

#Run ASRs
ASR_symbiont_type_AM_NMAM_lumped_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_AM_NMAM_lumped_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_AM_NMAM_lumped_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_symbiont_type_AM_NMAM_lumped_ER_yang")
load("./Output/ASR_symbiont_type_AM_NMAM_lumped_ARD_yang")
load("./Output/ASR_symbiont_type_AM_NMAM_lumped_SYM_yang")

save(ASR_symbiont_type_AM_NMAM_lumped_ER_yang, file = "./Output/ASR_symbiont_type_AM_NMAM_lumped_ER_yang")
save(ASR_symbiont_type_AM_NMAM_lumped_ARD_yang, file = "./Output/ASR_symbiont_type_AM_NMAM_lumped_ARD_yang")
save(ASR_symbiont_type_AM_NMAM_lumped_SYM_yang, file = "./Output/ASR_symbiont_type_AM_NMAM_lumped_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_AM_NMAM_lumped_ER_yang$AICc,
    ASR_symbiont_type_AM_NMAM_lumped_ARD_yang$AICc,
    ASR_symbiont_type_AM_NMAM_lumped_SYM_yang$AICc
  )
)

#ARD by far the best
ASR_symbiont_type_AM_NMAM_lumped_ARD_yang
plotMKmodel(ASR_symbiont_type_AM_NMAM_lumped_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_symbiont_type_AM_NMAM_lumped <-
  analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped %>%
  dplyr::select(Symbiotic_type)
row.names(dat_plot_symbiont_type_AM_NMAM_lumped) <-
  analysis_dat_CSR_symb_ASR_symbiont_type_AM_NMAM_lumped$Species_name
# dat_plot_symbiont_type$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_type$Symbiotic_type))
head(dat_plot_symbiont_type_AM_NMAM_lumped)

#Plot
pdf("./Output/ASRsymbiont_type_AM_NMAM_lumped.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_AM_NMAM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_AM_NMAM_lumped_ARD_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_AM_NMAM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_AM_NMAM_lumped_ARD_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()

####Second, lump all AM
#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped$Symbiotic_type<-
  ifelse(analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped$Symbiotic_type %in% c("AM","NMAM","AMNod","EcMAM"),
         "All_AM",analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped$Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped)
table(analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped$Symbiotic_type)

#Run ASRs
ASR_symbiont_type_All_AM_lumped_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_AM_lumped_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_AM_lumped_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_symbiont_type_All_AM_lumped_ER_yang")
load("./Output/ASR_symbiont_type_All_AM_lumped_ARD_yang")
load("./Output/ASR_symbiont_type_All_AM_lumped_SYM_yang")

save(ASR_symbiont_type_All_AM_lumped_ER_yang, file = "./Output/ASR_symbiont_type_All_AM_lumped_ER_yang")
save(ASR_symbiont_type_All_AM_lumped_ARD_yang, file = "./Output/ASR_symbiont_type_All_AM_lumped_ARD_yang")
save(ASR_symbiont_type_All_AM_lumped_SYM_yang, file = "./Output/ASR_symbiont_type_All_AM_lumped_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_All_AM_lumped_ER_yang$AICc,
    ASR_symbiont_type_All_AM_lumped_ARD_yang$AICc,
    ASR_symbiont_type_All_AM_lumped_SYM_yang$AICc
  )
)

#SYM by far the best
ASR_symbiont_type_All_AM_lumped_SYM_yang
plotMKmodel(ASR_symbiont_type_All_AM_lumped_SYM_yang)

#Plot them

#Create a data frame to plot the trait data
dat_plot_symbiont_type_All_AM_lumped <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped %>%
  dplyr::select(Symbiotic_type)
row.names(dat_plot_symbiont_type_All_AM_lumped) <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_AM_lumped$Species_name
# dat_plot_symbiont_type$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_type$Symbiotic_type))
head(dat_plot_symbiont_type_All_AM_lumped)

#Plot
pdf("./Output/ASRsymbiont_type_All_AM_lumped.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_AM_lumped_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_AM_lumped_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()


#####Third, lump all non-AM
#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped$Symbiotic_type<-
  ifelse(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped$Symbiotic_type %in% c("AM"),
         "Strict_AM","non_AM")
head(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped)
table(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped$Symbiotic_type)

#Run ASRs
ASR_symbiont_type_All_non_AM_lumped_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_non_AM_lumped_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_non_AM_lumped_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_symbiont_type_All_non_AM_lumped_ER_yang")
load("./Output/ASR_symbiont_type_All_non_AM_lumped_ARD_yang")
load("./Output/ASR_symbiont_type_All_non_AM_lumped_SYM_yang")

save(ASR_symbiont_type_All_non_AM_lumped_ER_yang, file = "./Output/ASR_symbiont_type_All_non_AM_lumped_ER_yang")
save(ASR_symbiont_type_All_non_AM_lumped_ARD_yang, file = "./Output/ASR_symbiont_type_All_non_AM_lumped_ARD_yang")
save(ASR_symbiont_type_All_non_AM_lumped_SYM_yang, file = "./Output/ASR_symbiont_type_All_non_AM_lumped_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_All_non_AM_lumped_ER_yang$AICc,
    ASR_symbiont_type_All_non_AM_lumped_ARD_yang$AICc,
    ASR_symbiont_type_All_non_AM_lumped_SYM_yang$AICc
  )
)

#ARD by far the best
ASR_symbiont_type_All_non_AM_lumped_ARD_yang
plotMKmodel(ASR_symbiont_type_All_non_AM_lumped_ARD_yang)

#Plot them

#Create a data frame to plot the trait data
dat_plot_symbiont_type_All_non_AM_lumped <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped %>%
  dplyr::select(Symbiotic_type)
row.names(dat_plot_symbiont_type_All_non_AM_lumped) <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_lumped$Species_name
# dat_plot_symbiont_type$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_type$Symbiotic_type))
head(dat_plot_symbiont_type_All_non_AM_lumped)
table(dat_plot_symbiont_type_All_non_AM_lumped$Symbiotic_type)

#Plot
pdf("./Output/ASRsymbiont_type_All_non_AM_lumped.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_non_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_non_AM_lumped_ARD_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_non_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_non_AM_lumped_ARD_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()

#####Fourth, lump all non-AM And all with AnyAM
#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped$Symbiotic_type<-
  ifelse(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped$Symbiotic_type %in% c("AM","NMAM","AMNod","EcMAM"),
         "Any_AM","non_AM")
head(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped)
table(analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped$Symbiotic_type)

#Run ASRs
ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ER_yang")
load("./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ARD_yang")
load("./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang")

save(ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ER_yang, file = "./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ER_yang")
save(ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ARD_yang, file = "./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ARD_yang")
save(ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang, file = "./Output/ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ER_yang$AICc,
    ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_ARD_yang$AICc,
    ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang$AICc
  )
)

#SYM or ER by far the best - they are the same in this case, because only two states. 
ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang
plotMKmodel(ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang)

#Plot them

#Create a data frame to plot the trait data
dat_plot_symbiont_type_All_non_AM_vs_any_AM_lumped <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped %>%
  dplyr::select(Symbiotic_type)
row.names(dat_plot_symbiont_type_All_non_AM_vs_any_AM_lumped) <-
  analysis_dat_CSR_symb_ASR_symbiont_type_All_non_AM_vs_any_AM_lumped$Species_name
# dat_plot_symbiont_type$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_type$Symbiotic_type))
head(dat_plot_symbiont_type_All_non_AM_vs_any_AM_lumped)

#Plot
pdf("./Output/ASRsymbiont_type_All_non_AM_vs_any_AM_lumped.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_non_AM_vs_any_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_type_All_non_AM_vs_any_AM_lumped,
  cols = list(Symbiotic_type = brewer.pal(n = 8, "Set2")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_symbiont_type_All_non_AM_vs_any_AM_lumped_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.3)
add.scale.bar()



############## CSR ASR

#Ways to treat CSR analytically:
#1. Turn it into a categorical variable: assign one of three categorical states based on what it is most selected for.
#2. Treat as three distinct continuous variables, and repeat the analyses for each.
#Neither of these captures perfectly what we are trying to measure, but they may do for our purposes. Let's have a look.
#I'll try only 1 for now, because it seems simplest and most informative.

analysis_dat_CSR_symb$selection_type <-
  apply(
    analysis_dat_CSR_symb %>% dplyr::select(C.selection, S.selection, R.selection),
    1,
    which.max
  )
analysis_dat_CSR_symb$selection_type <-
  gsub(
    pattern = "1",
    replacement = "C",
    x = analysis_dat_CSR_symb$selection_type
  )
analysis_dat_CSR_symb$selection_type <-
  gsub(
    pattern = "2",
    replacement = "S",
    x = analysis_dat_CSR_symb$selection_type
  )
analysis_dat_CSR_symb$selection_type <-
  gsub(
    pattern = "3",
    replacement = "R",
    x = analysis_dat_CSR_symb$selection_type
  )
table(analysis_dat_CSR_symb$selection_type) #Pretty equal numbers of all three types.

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_selection_type <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, selection_type)
head(analysis_dat_CSR_symb_ASR_selection_type)

#Run ASRs
ASR_selection_type_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_selection_type_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_selection_type_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_selection_type_ER_yang")
load("./Output/ASR_selection_type_ARD_yang")
load("./Output/ASR_selection_type_SYM_yang")

save(ASR_selection_type_ER_yang, file = "./Output/ASR_selection_type_ER_yang")
save(ASR_selection_type_ARD_yang, file = "./Output/ASR_selection_type_ARD_yang")
save(ASR_selection_type_SYM_yang, file = "./Output/ASR_selection_type_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_selection_type_ER_yang$AICc,
    ASR_selection_type_ARD_yang$AICc,
    ASR_selection_type_SYM_yang$AICc
  )
)

#ARD by far the best
ASR_selection_type_ARD_yang
plotMKmodel(ASR_selection_type_ARD_yang)
table(analysis_dat_CSR_symb$selection_type) #States are numbered in the modeling: this is what types the numbers represent, they are ordered aphabetically, it sems.

#Create a data frame to plot the trait data
dat_plot_selection_type <-
  analysis_dat_CSR_symb_ASR_selection_type %>%
  dplyr::select(selection_type)
row.names(dat_plot_selection_type) <-
  analysis_dat_CSR_symb_ASR_selection_type$Species_name
dat_plot_selection_type$selection_type <-
  as.numeric(as.factor(dat_plot_selection_type$selection_type))
head(dat_plot_selection_type)

#CSR ASR - Plot to Pdf
pdf("./Output/ASRCSRType.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_selection_type,
  cols = list(selection_type = brewer.pal(n = 3, "Accent")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_selection_type_ARD_yang$states,
           piecol = brewer.pal(n = 3, "Accent"),
           cex = 0.3)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_selection_type,
  cols = list(selection_type = brewer.pal(n = 3, "Accent")),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_selection_type_ARD_yang$states,
           piecol = brewer.pal(n = 3, "Accent"),
           cex = 0.3)
add.scale.bar()

#Ok, what do we see here? Reconstruction of the three CSR types, treated as three discrete categories.
#Some things that spring to mind:
#1. Looking at the distribution of the trait on the outside bar, there doesn't seem to be massive phylogenetic signal (we could quantify this)
#2. Or in other words, the three types are all mixed up, mostly.
#3  As a result throughout most of evolutionary history, it seems to be the case that we can't do mutch better than reconstruct 1/3, 1/3, 1/3 for the three types.
#4. Thus, evolution of distinct CSR-types is even more tippy than symbiont.
#5. Or, in other words. They generally don't precede, but follow symbiont switches.
#6. Alternative explantion CSR evolution happens at a much much faster rate. I.e. we just can't look back that far with ASR.
#7. Or, it's the method of splitting into 3 discrete categories, which creates the inevitbale simplification.
#8. Options: if we have fossile/other evidence of ancient CSR-types, we could fix some nodes and create a better ASR.
# Does this exist?


############## CSR ASR - Approach 2

table(analysis_dat_CSR_symb$CSR_categorical_level) #Pretty equal numbers of all three types.
#Turn it into a binary with CSR vs other
analysis_dat_CSR_symb$CSR_binary <-
  ifelse(
    analysis_dat_CSR_symb$CSR_categorical_level %in% c("C/CSR", "CR/CSR", "CS/CSR", "CSR", "R/CSR", "S/CSR", "SR/CSR"),
    "AnyCSR",
    "NoCSR"
  )
table(analysis_dat_CSR_symb$CSR_binary)

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_selection_type_binary <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, CSR_binary)
head(analysis_dat_CSR_symb_ASR_selection_type_binary)

#Run ASRs
ASR_selection_type_binary_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type_binary,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = c(1,0), #Fix root to any CSR - email Marco mid September 2020
    nstarts = 10,
    n.cores = 7
  )
ASR_selection_type_binary_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type_binary,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = c(1,0),
    nstarts = 10,
    n.cores = 7
  )
ASR_selection_type_binary_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_selection_type_binary,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = c(1,0),
    nstarts = 10,
    n.cores = 7
  )

#Save all model ran
load("./Output/ASR_selection_type_binary_ER_yang")
load("./Output/ASR_selection_type_binary_ARD_yang")
load("./Output/ASR_selection_type_binary_SYM_yang")

save(ASR_selection_type_binary_ER_yang, file = "./Output/ASR_selection_type_binary_ER_yang")
save(ASR_selection_type_binary_ARD_yang, file = "./Output/ASR_selection_type_binary_ARD_yang")
save(ASR_selection_type_binary_SYM_yang, file = "./Output/ASR_selection_type_binary_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_selection_type_binary_ER_yang$AICc,
    ASR_selection_type_binary_ARD_yang$AICc,
    ASR_selection_type_binary_SYM_yang$AICc
  )
)

#ARD by far the best
ASR_selection_type_binary_ARD_yang
plotMKmodel(ASR_selection_type_binary_ARD_yang)
table(analysis_dat_CSR_symb_ASR_selection_type_binary$CSR_binary) #States are numbered in the modeling: this is what types the numbers represent, they are ordered aphabetically, it sems.

# #Create a data frame to plot the trait data
dat_plot_selection_type_binary <-
   analysis_dat_CSR_symb_ASR_selection_type_binary %>%
   dplyr::select(CSR_binary)
 row.names(dat_plot_selection_type_binary) <-
   analysis_dat_CSR_symb_ASR_selection_type_binary$Species_name
  dat_plot_selection_type_binary$CSR_binary <-
    as.numeric(as.factor(dat_plot_selection_type_binary$CSR_binary))
 head(dat_plot_selection_type_binary)

 #CSR ASR - Plot to Pdf
pdf("./Output/ASRCSRType_binary.pdf",
        width = 20,
        height = 20)
    trait.plot(
      tree = analysis_tree,
      dat = dat_plot_selection_type_binary,
      cols = list(CSR_binary = brewer.pal(n = 3, "Accent")),
      type = "f",
      legend = T,
      w = 1 / 40,
      edge.width = 2,
      cex.lab = 0.01,
      tip.color = "white",
      show.node.label = T
    )
    nodelabels(pie = ASR_selection_type_binary_ARD_yang$states,
               piecol = brewer.pal(n = 3, "Accent"),
               cex = 0.3)
    add.scale.bar()
    dev.off()

    #Plot to screen
    trait.plot(
      tree = analysis_tree,
      dat = dat_plot_selection_type_binary,
      cols = list(CSR_binary = brewer.pal(n = 3, "Accent")),
      type = "f",
      legend = T,
      w = 1 / 40,
      edge.width = 2,
      cex.lab = 0.01,
      tip.color = "white",
      show.node.label = T
    )
    nodelabels(pie = ASR_selection_type_binary_ARD_yang$states,
               piecol = brewer.pal(n = 3, "Accent"),
               cex = 0.3)
    add.scale.bar()

######Correlated evolution between the two variables

#Let's run a combined model, modelling to traits simultaneously first.

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type, CSR_binary)
head(analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary)

#Run ASRs
ASR_symbiont_selection_type_binary_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_selection_type_binary_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_symbiont_selection_type_binary_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary,
    rate.cat = 1,
    model = "SYM",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )

#

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_selection_type_binary_ER_yang$AICc,
    ASR_symbiont_selection_type_binary_ARD_yang$AICc,
    ASR_symbiont_selection_type_binary_SYM_yang$AICc
  )
)



#Save all model ran
load("./Output/ASR_symbiont_selection_type_binary_ER_yang")
load("./Output/ASR_symbiont_selection_type_binary_ARD_yang")
load("./Output/ASR_symbiont_selection_type_binary_SYM_yang")

save(ASR_symbiont_selection_type_binary_ER_yang, file = "./Output/ASR_symbiont_selection_type_binary_ER_yang")
save(ASR_symbiont_selection_type_binary_ARD_yang, file = "./Output/ASR_symbiont_selection_type_binary_ARD_yang")
save(ASR_symbiont_selection_type_binary_SYM_yang, file = "./Output/ASR_symbiont_selection_type_binary_SYM_yang")

#Let's look at the best (=ARD) model
ASR_symbiont_selection_type_binary_ARD_yang
plotMKmodel(ASR_symbiont_selection_type_binary_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_symbiont_selection_type_binary <-
  analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary %>%
  dplyr::select(Symbiotic_type, CSR_binary)
row.names(dat_plot_symbiont_selection_type_binary) <-
  analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$Species_name
# dat_plot_symbiont_selection_type_binary$Symbiotic_type <-
#   as.numeric(as.factor(dat_plot_symbiont_selection_type_binary$Symbiotic_type))
# dat_plot_symbiont_selection_type_binary$CSR_binary <-
#   as.numeric(as.factor(dat_plot_symbiont_selection_type_binary$CSR_binary))
head(dat_plot_symbiont_selection_type_binary)

#CSR ASR - plot pdf
pdf("./Output/ASRSymbiontCSRTypeBinary.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type_binary,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    CSR_binary = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(
  pie = ASR_symbiont_selection_type_binary_ER_yang$states,
  piecol = c(brewer.pal(n = 11, "Set3"), brewer.pal(n = 11, "Paired")),
  cex = 0.3
)
legend(
  legend = paste(
    rep(
      unique(
        analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$Symbiotic_type
      ),
      each = length(
        unique(
          analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$CSR_binary
        )
      )
    ),
    unique(
      analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$CSR_binary
    ),
    sep = " & "
  ),
  x = "bottomright",
  fill = c(brewer.pal(n = 11, "Set3"), brewer.pal(n = 11, "Paired"))
)
add.scale.bar()
dev.off()

trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type_binary,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    CSR_binary = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(
  pie = ASR_symbiont_selection_type_binary_ER_yang$states,
  piecol = c(brewer.pal(n = 11, "Set3"), brewer.pal(n = 11, "Paired")),
  cex = 0.3
)
legend(
  legend = paste(
    rep(
      unique(
        analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$Symbiotic_type
      ),
      each = length(
        unique(
          analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$CSR_binary
        )
      )
    ),
    unique(
      analysis_dat_CSR_symb_ASR_symbiont_selection_type_binary$CSR_binary
    ),
    sep = " & "
  ),
  x = "bottomright",
  fill = c(brewer.pal(n = 11, "Set3"), brewer.pal(n = 11, "Paired"))
)
add.scale.bar()

#Ok, what are we seeing here?
#A joint reconstruction of both categorical variables simultaneously.
#This gets at your quesion most directly.
#What does it tell us? Well something a bit more complicated than before, and also partally contradictory, I would say.
#For instance from the ancestral state of angiosperms AM with S-type, we first get a transition towards AM & R
#Or, this is CSR-selection shift first, then symbiont.
#But, there's lots and lots of detail in this figure, so this may not be true universally.
#If you want to present something like this in the paper, we would need to think of a clever colouring scheme.
#You could work with hues of given colour to indicate CSR-shifts within a symbiont type (or the other way around) - could be quite beautiful, I think.
#Big caveat: This model is not reliable yet. Only the simples version (Equal Rate - ER) has finished so far.
#For a complex evolutionary scenario like this, these tend to not be great.
#My idea expanding this line of analysis is probably fruitful.
#Things to decide to that:
#1. Do we manually update species? or accept the 90% overlap (my thought)
#2. Is the categorical approach for CSR acceptable?
#3. Do we have justified grounds to fix some nodes, particularly for the CSR-variable?


###############Quantitative analysis for CSR

#Reconstruct the CSR spectrum as three separate quantitative variables

#Data formatting create the vectors
vector_C_selection <- analysis_dat_CSR_symb$C.selection
names(vector_C_selection) <- analysis_dat_CSR_symb$Species_name
vector_S_selection <- analysis_dat_CSR_symb$S.selection
names(vector_S_selection) <- analysis_dat_CSR_symb$Species_name
vector_R_selection <- analysis_dat_CSR_symb$R.selection
names(vector_R_selection) <- analysis_dat_CSR_symb$Species_name
#For ease of plotting, order vector same order as in tree
vector_C_selection <-
  vector_C_selection[match(analysis_tree$tip.label, names(vector_C_selection))]
vector_S_selection <-
  vector_S_selection[match(analysis_tree$tip.label, names(vector_S_selection))]
vector_R_selection <-
  vector_R_selection[match(analysis_tree$tip.label, names(vector_R_selection))]

head(vector_C_selection)
head(vector_S_selection)
head(vector_R_selection)
#Look all good.

#Run the models
quant_ASR_C_selection <-
  anc.recon(trait_data = vector_C_selection, tree = analysis_tree)
save(quant_ASR_C_selection, file = "./Output/quant_ASR_C_selection")

quant_ASR_S_selection <-
  anc.recon(trait_data = vector_S_selection, tree = analysis_tree)
save(quant_ASR_S_selection, file = "./Output/quant_ASR_S_selection")

quant_ASR_R_selection <-
  anc.recon(trait_data = vector_R_selection, tree = analysis_tree)
save(quant_ASR_R_selection, file = "./Output/quant_ASR_R_selection")

head(quant_ASR_C_selection)
head(quant_ASR_S_selection)
head(quant_ASR_R_selection)

#Plot the three reconstruction, and overlay the symbiotic reconstructions on top.

#Symbionts ASR with quantitative ASR of CSelection - plot to pdf
pdf(
  "./Output/ASRSymbiontType_QuantCSelection.pdf",
  width = 20,
  height = 20
)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_C_selection[match(analysis_tree$edge[, 1], names(quant_ASR_C_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "C Selection",
  prompt = F,
  lims = c(min(quant_ASR_C_selection), max(quant_ASR_C_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_C_selection[match(analysis_tree$edge[, 1], names(quant_ASR_C_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "C Selection",
  prompt = F,
  lims = c(min(quant_ASR_C_selection), max(quant_ASR_C_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()


#Symbionts ASR with quantitative ASR of S_Selection - plot to pdf
pdf(
  "./Output/ASRSymbiontType_QuantSSelection.pdf",
  width = 20,
  height = 20
)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_S_selection[match(analysis_tree$edge[, 1], names(quant_ASR_S_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "S Selection",
  prompt = F,
  lims = c(min(quant_ASR_S_selection), max(quant_ASR_S_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_S_selection[match(analysis_tree$edge[, 1], names(quant_ASR_S_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "S Selection",
  prompt = F,
  lims = c(min(quant_ASR_S_selection), max(quant_ASR_S_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()


#Symbionts ASR with quantitative ASR of R_Selection - plot to pdf
pdf(
  "./Output/ASRSymbiontType_QuantRSelection.pdf",
  width = 20,
  height = 20
)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_R_selection[match(analysis_tree$edge[, 1], names(quant_ASR_R_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "R Selection",
  prompt = F,
  lims = c(min(quant_ASR_R_selection), max(quant_ASR_R_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()
dev.off()

#Plot to screen
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_symbiont_selection_type,
  cols = list(
    Symbiotic_type = brewer.pal(n = 8, "Set2"),
    selection_type = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = T,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T,
  edge.color = inferno(100)[cut(quant_ASR_R_selection[match(analysis_tree$edge[, 1], names(quant_ASR_R_selection[, 1])), 1], breaks =
                                  100)]
)
nodelabels(pie = ASR_symbiont_type_SYM_yang$states,
           piecol = brewer.pal(n = 8, "Set2"),
           cex = 0.15)
add.color.bar(
  100,
  inferno(100),
  title = "R Selection",
  prompt = F,
  lims = c(min(quant_ASR_R_selection), max(quant_ASR_R_selection)),
  fsize = 0.8,
  x = -100,
  y = -50
)
add.scale.bar()

#Ok so what are we seeing here.
#It's the categorical recosntructions of symbionts plotted onto the coloured tree branches, indicating C/S/R selection
#General picture it seems: CSR-levels are mostly moderate throughout evolution, and big shifts are 'tippy'.
#That means either (1) biological result, the distinct strategies evolved only late in evolution.
#Or (2) it means we are dealing with an artifact.
#What I mean with 2: perhaps we just can't really evaluate the ancient evolution of CSR based on only tip data, because it simply evolves to fast.
#here too, the solution would be fixing nodes based on fossile evidence, I guess?






###Potential follow-up
#ASRs of the three CSR values. 3x
#Plot these on top of each other
#Pagel's model, need two binary variables.
