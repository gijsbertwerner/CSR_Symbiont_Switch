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

smith_brown_tree <- read.tree("./Data/ALLMB.tre")
smith_brown_tree

# Data Cleaning -----------------------------------------------------------

####Clean data file

#Formatting of levels
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

#Species formatting
head(dat_CSR_symb$Species_name)
head(smith_brown_tree$tip.label)
dat_CSR_symb$Species_name <- gsub(pattern = " ",
                                  replacement = "_",
                                  dat_CSR_symb$Species_name)

#How many of the 3014 species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb$Species_name, smith_brown_tree$tip.label))
#Only 320 of 3014 lacking. I.e. ~90% is present. That seems good enough.

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
#Marco manually checked the above two files and suggested name substitutions (wrong synonyms etc.)
###Substitutions of data species - based on manual checking by Marco, happens here.
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

#### Clean phylogeny

#Extract the appropriate subtree from Smith&Brown
analysis_tree <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb$Species_name)
  )
analysis_tree

#PLot big analysis tree for finding species
pdf("./Output/BigAnalysisTree.pdf",
    width = 30,
    height = 30)
plot.phylo(
  x = analysis_tree,
  type = "f",
  show.tip.label = T,
  cex = 0.1,
  label.offset = 0.2
)
dev.off()

#### Match dataset and phylogeny

analysis_dat_CSR_symb <-
  dat_CSR_symb %>% filter(Species_name %in% analysis_tree$tip.label)
nrow(analysis_dat_CSR_symb)
#More species in data set, than in tree.
#This should not be possible. Could be due to duplicates in the dataset?
analysis_dat_CSR_symb[duplicated(analysis_dat_CSR_symb$Species_name),]
#Yes, some duplicates, but only about ~30. Small variations in CSR-values, agreement on symbiont type

#Simply remove the duplicates.
analysis_dat_CSR_symb <-
  analysis_dat_CSR_symb[!duplicated(analysis_dat_CSR_symb$Species_name), ]
nrow(analysis_dat_CSR_symb)
#Now matches

# Defining variables for analyses -------------------------------------------------------

#Create the categorical variables to analyse

#####First for the symbionts
head(analysis_dat_CSR_symb)
table(analysis_dat_CSR_symb$Symbiotic_type)

#Turn the categorical variable in two binary ones.
analysis_dat_CSR_symb$Symbiotic_type_AMOnly_Rest <-
  ifelse(analysis_dat_CSR_symb$Symbiotic_type %in% c("AM"),
         "AMOnly",
         "Rest")
table(analysis_dat_CSR_symb$Symbiotic_type_AMOnly_Rest)
analysis_dat_CSR_symb$Symbiotic_type_AnyAM_Rest <-
  ifelse(
    analysis_dat_CSR_symb$Symbiotic_type %in% c("AM", "NMAM", "AMNod", "EcMAM"),
    "AnyAM",
    "Rest"
  )
table(analysis_dat_CSR_symb$Symbiotic_type_AnyAM_Rest)

######And second for the CSR strategies
table(analysis_dat_CSR_symb$CSR_categorical_level) #Pretty equal numbers of all three types.

#Create a number of alternative CSR-cutoffs, using cosine similarity (email Marco October 30, 2020).
analysis_dat_CSR_symb$CSR_binary_70 <-
  ifelse(
    analysis_dat_CSR_symb$CSR_categorical_level %in% c(
      "C/CSR",
      "CR/CSR",
      "CS/CSR",
      "CSR",
      "R/CSR",
      "S/CSR",
      "SR/CSR",
      "CR",
      "CS",
      "SR",
      "C/CR",
      "C/CS",
      "R/CR",
      "R/SR",
      "S/CS",
      "S/SR"
    ),
    "CSR70",
    "NoCSR"
  )
analysis_dat_CSR_symb$CSR_binary_85 <-
  ifelse(
    analysis_dat_CSR_symb$CSR_categorical_level %in% c(
      "C/CSR",
      "CR/CSR",
      "CS/CSR",
      "CSR",
      "R/CSR",
      "S/CSR",
      "SR/CSR",
      "CR",
      "CS",
      "SR"
    ),
    "CSR85",
    "NoCSR"
  )
analysis_dat_CSR_symb$CSR_binary_90 <-
  ifelse(
    analysis_dat_CSR_symb$CSR_categorical_level %in% c("C/CSR", "CR/CSR", "CS/CSR", "CSR", "R/CSR", "S/CSR", "SR/CSR"),
    "CSR90",
    "NoCSR"
  )
analysis_dat_CSR_symb$CSR_binary_92 <-
  ifelse(
    analysis_dat_CSR_symb$CSR_categorical_level %in% c("CR/CSR", "CS/CSR", "CSR", "SR/CSR"),
    "CSR92",
    "NoCSR"
  )
analysis_dat_CSR_symb$CSR_binary_95 <-
  ifelse(analysis_dat_CSR_symb$CSR_categorical_level %in% c("CSR"),
         "CSR95",
         "NoCSR")
table(analysis_dat_CSR_symb$CSR_binary_70)
table(analysis_dat_CSR_symb$CSR_binary_85)
table(analysis_dat_CSR_symb$CSR_binary_90)
table(analysis_dat_CSR_symb$CSR_binary_92)
table(analysis_dat_CSR_symb$CSR_binary_95)

##Check data file
head(analysis_dat_CSR_symb)

write.csv(
  analysis_dat_CSR_symb,
  quote = F,
  row.names = F,
  file = "./Output/AnalysedData.csv"
)


# Analyses ----------------------------------------------------------------


###### Analysing the symbiotic states

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type)
table(analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type<-as.numeric(as.factor(analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type))
table(analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type)

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
save(ASR_symbiont_type_ER_yang, file = "./Output/ASR_symbiont_type_ER_yang")
save(ASR_symbiont_type_ARD_yang, file = "./Output/ASR_symbiont_type_ARD_yang")
save(ASR_symbiont_type_SYM_yang, file = "./Output/ASR_symbiont_type_SYM_yang")

load("./Output/ASR_symbiont_type_ER_yang")
load("./Output/ASR_symbiont_type_ARD_yang")
load("./Output/ASR_symbiont_type_SYM_yang")

#Which is the best model, using AIC-criteria?
akaike.weights(
  c(
    ASR_symbiont_type_ER_yang$AICc,
    ASR_symbiont_type_ARD_yang$AICc,
    ASR_symbiont_type_SYM_yang$AICc
  )
)

#ARD is the best
ASR_symbiont_type_SYM_yang
plotMKmodel(ASR_symbiont_type_SYM_yang)
table(analysis_dat_CSR_symb$Symbiotic_type) #States are numbered in the modeling: this is what types the numbers represent, they are ordered aphabetically, it sems.

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
legend(legend=names(table(analysis_dat_CSR_symb$Symbiotic_type)),
       x = "bottomright",
       fill = brewer.pal(n = 8, "Set2"),
       cex = 2)
add.scale.bar()
dev.off()


######Correlated evolution between the two variables

######90

#First we'll study the baseline case (cutoff cosine similarity = 90, AnyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AnyAM_Rest, CSR_binary_90)
head(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90
)
states_AnyAM_Rest_CSR_binary_90 <-
  c("AnyAM & CSR90", "AnyAM & NoCSR", "noAM & CSR90", "noAM & NoCSR")

analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type_AnyAM_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type_AnyAM_Rest
    )
  )
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90
  ))
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$CSR_binary_90
)

#Run ASRs
ASR_AnyAM_Rest_CSR_binary_90_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_90_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_90_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90,
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
    ASR_AnyAM_Rest_CSR_binary_90_ER_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_90_ARD_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_90_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AnyAM_Rest_CSR_binary_90_ER_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_90_ER_yang")
save(ASR_AnyAM_Rest_CSR_binary_90_ARD_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_90_ARD_yang")
save(ASR_AnyAM_Rest_CSR_binary_90_SYM_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_90_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AnyAM_Rest_CSR_binary_90_ARD_yang
plotMKmodel(ASR_AnyAM_Rest_CSR_binary_90_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AnyAM_Rest_CSR_binary_90 <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90 %>%
  dplyr::select(Symbiotic_type_AnyAM_Rest, CSR_binary_90)
row.names(dat_plot_AnyAM_Rest_CSR_binary_90) <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_90$Species_name
head(dat_plot_AnyAM_Rest_CSR_binary_90)

plotvec_symbiont_binary_selection_type_binary<-
  c("#e31a1c","#fb9a99",     #Red is AMF
    "#8c510a","#d8b365",     #Brown is AMNod
    "#1f78b4","#a6cee3",      #Blue is ECM
    "#a6cee3","#cab2d6",     #Purple is ECMAM
    "#33a02c","#b2df8a",      #Green is ERM  
    "#fec44f","#fee391",     #Yellow is NM
    "#525252","#bdbdbd",    #Grey is NMAM
    "#014636" ,"#02818a")    #Turqouise is OM 

#CSR ASR - plot pdf
pdf("./Output/ASR_AnyAM_Rest_CSR_binary_90.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AnyAM_Rest_CSR_binary_90,
  cols = list(
    Symbiotic_type_AnyAM_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_90 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AnyAM_Rest_CSR_binary_90_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AnyAM_Rest_CSR_binary_90,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)","Rest (EcM,ErM,NM,OM)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_90$Symbiotic_type_AnyAM_Rest
vec_selection_type_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_90$CSR_binary_90
names(vec_symbiont_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_90)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_90)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AnyAM_Rest_CSR_binary_90_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AnyAM_Rest_CSR_binary_90_ARD
plot(pagel_AnyAM_Rest_CSR_binary_90_ARD)

pdf("./Output/Plot_pagel_AnyAM_Rest_CSR_binary_90_ARD.pdf")
plot(pagel_AnyAM_Rest_CSR_binary_90_ARD)
dev.off()

save.image()

#Second, the baseline case with OnlyAM (cutoff cosine similarity = 90, OnlyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AMOnly_Rest, CSR_binary_90)
head(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90
)
states_AMOnly_Rest_CSR_binary_90 <-
  c("AMOnly & CSR90", "AMOnly & NoCSR", "noAM & CSR90", "noAM & NoCSR")

analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type_AMOnly_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type_AMOnly_Rest
    )
  )
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90
  ))
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$CSR_binary_90
)

#Run ASRs
ASR_AMOnly_Rest_CSR_binary_90_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_90_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_90_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90,
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
    ASR_AMOnly_Rest_CSR_binary_90_ER_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_90_ARD_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_90_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AMOnly_Rest_CSR_binary_90_ER_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_90_ER_yang")
save(ASR_AMOnly_Rest_CSR_binary_90_ARD_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_90_ARD_yang")
save(ASR_AMOnly_Rest_CSR_binary_90_SYM_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_90_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AMOnly_Rest_CSR_binary_90_ARD_yang
plotMKmodel(ASR_AMOnly_Rest_CSR_binary_90_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AMOnly_Rest_CSR_binary_90 <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90 %>%
  dplyr::select(Symbiotic_type_AMOnly_Rest, CSR_binary_90)
row.names(dat_plot_AMOnly_Rest_CSR_binary_90) <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_90$Species_name
head(dat_plot_AMOnly_Rest_CSR_binary_90)

#CSR ASR - plot pdf
pdf("./Output/ASR_AMOnly_Rest_CSR_binary_90.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AMOnly_Rest_CSR_binary_90,
  cols = list(
    Symbiotic_type_AMOnly_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_90 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AMOnly_Rest_CSR_binary_90_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AMOnly_Rest_CSR_binary_90,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AMOnly (AM)","Rest (All other types)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_90$Symbiotic_type_AMOnly_Rest
vec_selection_type_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_90$CSR_binary_90
names(vec_symbiont_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_90)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_90)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AMOnly_Rest_CSR_binary_90_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AMOnly_Rest_CSR_binary_90_ARD
plot(pagel_AMOnly_Rest_CSR_binary_90_ARD)

pdf("./Output/Plot_pagel_AMOnly_Rest_CSR_binary_90_ARD.pdf")
plot(pagel_AMOnly_Rest_CSR_binary_90_ARD)
dev.off()

save.image()


# Sensitivity Analyses ----------------------------------------------------

######70

#First we'll study the baseline case (cutoff cosine similarity = 70, AnyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AnyAM_Rest, CSR_binary_70)
head(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70
)
states_AnyAM_Rest_CSR_binary_70 <-
  c("AnyAM & CSR70", "AnyAM & NoCSR", "noAM & CSR70", "noAM & NoCSR")

analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type_AnyAM_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type_AnyAM_Rest
    )
  )
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70
  ))
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$CSR_binary_70
)

#Run ASRs
ASR_AnyAM_Rest_CSR_binary_70_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_70_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_70_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70,
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
    ASR_AnyAM_Rest_CSR_binary_70_ER_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_70_ARD_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_70_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AnyAM_Rest_CSR_binary_70_ER_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_70_ER_yang")
save(ASR_AnyAM_Rest_CSR_binary_70_ARD_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_70_ARD_yang")
save(ASR_AnyAM_Rest_CSR_binary_70_SYM_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_70_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AnyAM_Rest_CSR_binary_70_ARD_yang
plotMKmodel(ASR_AnyAM_Rest_CSR_binary_70_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AnyAM_Rest_CSR_binary_70 <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70 %>%
  dplyr::select(Symbiotic_type_AnyAM_Rest, CSR_binary_70)
row.names(dat_plot_AnyAM_Rest_CSR_binary_70) <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_70$Species_name
head(dat_plot_AnyAM_Rest_CSR_binary_70)

#CSR ASR - plot pdf
pdf("./Output/ASR_AnyAM_Rest_CSR_binary_70.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AnyAM_Rest_CSR_binary_70,
  cols = list(
    Symbiotic_type_AnyAM_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_70 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AnyAM_Rest_CSR_binary_70_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AnyAM_Rest_CSR_binary_70,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)","Rest (EcM,ErM,NM,OM)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_70$Symbiotic_type_AnyAM_Rest
vec_selection_type_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_70$CSR_binary_70
names(vec_symbiont_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_70)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_70)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AnyAM_Rest_CSR_binary_70_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AnyAM_Rest_CSR_binary_70_ARD
plot(pagel_AnyAM_Rest_CSR_binary_70_ARD)

pdf("./Output/Plot_pagel_AnyAM_Rest_CSR_binary_70_ARD.pdf")
plot(pagel_AnyAM_Rest_CSR_binary_70_ARD)
dev.off()

save.image()

#Second, the baseline case with OnlyAM (cutoff cosine similarity = 70, OnlyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AMOnly_Rest, CSR_binary_70)
head(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70
)
states_AMOnly_Rest_CSR_binary_70 <-
  c("AMOnly & CSR70", "AMOnly & NoCSR", "noAM & CSR70", "noAM & NoCSR")

analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type_AMOnly_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type_AMOnly_Rest
    )
  )
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70
  ))
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$CSR_binary_70
)

#Run ASRs
ASR_AMOnly_Rest_CSR_binary_70_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_70_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_70_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70,
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
    ASR_AMOnly_Rest_CSR_binary_70_ER_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_70_ARD_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_70_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AMOnly_Rest_CSR_binary_70_ER_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_70_ER_yang")
save(ASR_AMOnly_Rest_CSR_binary_70_ARD_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_70_ARD_yang")
save(ASR_AMOnly_Rest_CSR_binary_70_SYM_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_70_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AMOnly_Rest_CSR_binary_70_ARD_yang
plotMKmodel(ASR_AMOnly_Rest_CSR_binary_70_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AMOnly_Rest_CSR_binary_70 <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70 %>%
  dplyr::select(Symbiotic_type_AMOnly_Rest, CSR_binary_70)
row.names(dat_plot_AMOnly_Rest_CSR_binary_70) <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_70$Species_name
head(dat_plot_AMOnly_Rest_CSR_binary_70)


#CSR ASR - plot pdf
pdf("./Output/ASR_AMOnly_Rest_CSR_binary_70.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AMOnly_Rest_CSR_binary_70,
  cols = list(
    Symbiotic_type_AMOnly_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_70 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AMOnly_Rest_CSR_binary_70_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AMOnly_Rest_CSR_binary_70,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AMOnly (AM)","Rest (All other types)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_70$Symbiotic_type_AMOnly_Rest
vec_selection_type_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_70$CSR_binary_70
names(vec_symbiont_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_70)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_70)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AMOnly_Rest_CSR_binary_70_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AMOnly_Rest_CSR_binary_70_ARD
plot(pagel_AMOnly_Rest_CSR_binary_70_ARD)

pdf("./Output/Plot_pagel_AMOnly_Rest_CSR_binary_70_ARD.pdf")
plot(pagel_AMOnly_Rest_CSR_binary_70_ARD)
dev.off()

save.image()

######85

#First we'll study the baseline case (cutoff cosine similarity = 85, AnyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AnyAM_Rest, CSR_binary_85)
head(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85
)
states_AnyAM_Rest_CSR_binary_85 <-
  c("AnyAM & CSR85", "AnyAM & NoCSR", "noAM & CSR85", "noAM & NoCSR")

analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type_AnyAM_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type_AnyAM_Rest
    )
  )
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85
  ))
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$CSR_binary_85
)

#Run ASRs
ASR_AnyAM_Rest_CSR_binary_85_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_85_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_85_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85,
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
    ASR_AnyAM_Rest_CSR_binary_85_ER_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_85_ARD_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_85_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AnyAM_Rest_CSR_binary_85_ER_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_85_ER_yang")
save(ASR_AnyAM_Rest_CSR_binary_85_ARD_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_85_ARD_yang")
save(ASR_AnyAM_Rest_CSR_binary_85_SYM_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_85_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AnyAM_Rest_CSR_binary_85_ARD_yang
plotMKmodel(ASR_AnyAM_Rest_CSR_binary_85_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AnyAM_Rest_CSR_binary_85 <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85 %>%
  dplyr::select(Symbiotic_type_AnyAM_Rest, CSR_binary_85)
row.names(dat_plot_AnyAM_Rest_CSR_binary_85) <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_85$Species_name
head(dat_plot_AnyAM_Rest_CSR_binary_85)

#CSR ASR - plot pdf
pdf("./Output/ASR_AnyAM_Rest_CSR_binary_85.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AnyAM_Rest_CSR_binary_85,
  cols = list(
    Symbiotic_type_AnyAM_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_85 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AnyAM_Rest_CSR_binary_85_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AnyAM_Rest_CSR_binary_85,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)","Rest (EcM,ErM,NM,OM)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_85$Symbiotic_type_AnyAM_Rest
vec_selection_type_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_85$CSR_binary_85
names(vec_symbiont_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_85)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_85)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AnyAM_Rest_CSR_binary_85_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AnyAM_Rest_CSR_binary_85_ARD
plot(pagel_AnyAM_Rest_CSR_binary_85_ARD)

pdf("./Output/Plot_pagel_AnyAM_Rest_CSR_binary_85_ARD.pdf")
plot(pagel_AnyAM_Rest_CSR_binary_85_ARD)
dev.off()

save.image()

#Second, the baseline case with OnlyAM (cutoff cosine similarity = 85, OnlyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AMOnly_Rest, CSR_binary_85)
head(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85
)
states_AMOnly_Rest_CSR_binary_85 <-
  c("AMOnly & CSR85", "AMOnly & NoCSR", "noAM & CSR85", "noAM & NoCSR")

analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type_AMOnly_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type_AMOnly_Rest
    )
  )
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85
  ))
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$CSR_binary_85
)

#Run ASRs
ASR_AMOnly_Rest_CSR_binary_85_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_85_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_85_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85,
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
    ASR_AMOnly_Rest_CSR_binary_85_ER_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_85_ARD_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_85_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AMOnly_Rest_CSR_binary_85_ER_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_85_ER_yang")
save(ASR_AMOnly_Rest_CSR_binary_85_ARD_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_85_ARD_yang")
save(ASR_AMOnly_Rest_CSR_binary_85_SYM_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_85_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AMOnly_Rest_CSR_binary_85_ARD_yang
plotMKmodel(ASR_AMOnly_Rest_CSR_binary_85_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AMOnly_Rest_CSR_binary_85 <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85 %>%
  dplyr::select(Symbiotic_type_AMOnly_Rest, CSR_binary_85)
row.names(dat_plot_AMOnly_Rest_CSR_binary_85) <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_85$Species_name
head(dat_plot_AMOnly_Rest_CSR_binary_85)

plotvec_symbiont_selection_type_binary_70<-
  c("#e31a1c","#fb9a99",     #Red is AMF
    "#8c510a","#d8b365",     #Brown is AMNod
    "#1f78b4","#a6cee3",      #Blue is ECM
    "#a6cee3","#cab2d6",     #Purple is ECMAM
    "#33a02c","#b2df8a",      #Green is ERM  
    "#fec44f","#fee391",     #Yellow is NM
    "#525252","#bdbdbd",    #Grey is NMAM
    "#014636" ,"#02818a")    #Turqouise is OM 

#CSR ASR - plot pdf
pdf("./Output/ASR_AMOnly_Rest_CSR_binary_85.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AMOnly_Rest_CSR_binary_85,
  cols = list(
    Symbiotic_type_AMOnly_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_85 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AMOnly_Rest_CSR_binary_85_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AMOnly_Rest_CSR_binary_85,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AMOnly (AM)","Rest (All other types)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_85$Symbiotic_type_AMOnly_Rest
vec_selection_type_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_85$CSR_binary_85
names(vec_symbiont_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_85)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_85)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AMOnly_Rest_CSR_binary_85_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AMOnly_Rest_CSR_binary_85_ARD
plot(pagel_AMOnly_Rest_CSR_binary_85_ARD)

pdf("./Output/Plot_pagel_AMOnly_Rest_CSR_binary_85_ARD.pdf")
plot(pagel_AMOnly_Rest_CSR_binary_85_ARD)
dev.off()

save.image()

######92

#First we'll study the baseline case (cutoff cosine similarity = 92, AnyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AnyAM_Rest, CSR_binary_92)
head(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92
)
states_AnyAM_Rest_CSR_binary_92 <-
  c("AnyAM & CSR92", "AnyAM & NoCSR", "noAM & CSR92", "noAM & NoCSR")

analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type_AnyAM_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type_AnyAM_Rest
    )
  )
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92
  ))
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$CSR_binary_92
)

#Run ASRs
ASR_AnyAM_Rest_CSR_binary_92_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_92_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_92_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92,
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
    ASR_AnyAM_Rest_CSR_binary_92_ER_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_92_ARD_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_92_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AnyAM_Rest_CSR_binary_92_ER_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_92_ER_yang")
save(ASR_AnyAM_Rest_CSR_binary_92_ARD_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_92_ARD_yang")
save(ASR_AnyAM_Rest_CSR_binary_92_SYM_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_92_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AnyAM_Rest_CSR_binary_92_ARD_yang
plotMKmodel(ASR_AnyAM_Rest_CSR_binary_92_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AnyAM_Rest_CSR_binary_92 <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92 %>%
  dplyr::select(Symbiotic_type_AnyAM_Rest, CSR_binary_92)
row.names(dat_plot_AnyAM_Rest_CSR_binary_92) <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_92$Species_name
head(dat_plot_AnyAM_Rest_CSR_binary_92)

#CSR ASR - plot pdf
pdf("./Output/ASR_AnyAM_Rest_CSR_binary_92.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AnyAM_Rest_CSR_binary_92,
  cols = list(
    Symbiotic_type_AnyAM_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_92 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AnyAM_Rest_CSR_binary_92_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AnyAM_Rest_CSR_binary_92,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)","Rest (EcM,ErM,NM,OM)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_92$Symbiotic_type_AnyAM_Rest
vec_selection_type_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_92$CSR_binary_92
names(vec_symbiont_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_92)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_92)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AnyAM_Rest_CSR_binary_92_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AnyAM_Rest_CSR_binary_92_ARD
plot(pagel_AnyAM_Rest_CSR_binary_92_ARD)

pdf("./Output/Plot_pagel_AnyAM_Rest_CSR_binary_92_ARD.pdf")
plot(pagel_AnyAM_Rest_CSR_binary_92_ARD)
dev.off()

save.image()

#Second, the baseline case with OnlyAM (cutoff cosine similarity = 92, OnlyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AMOnly_Rest, CSR_binary_92)
head(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92
)
states_AMOnly_Rest_CSR_binary_92 <-
  c("AMOnly & CSR92", "AMOnly & NoCSR", "noAM & CSR92", "noAM & NoCSR")

analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type_AMOnly_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type_AMOnly_Rest
    )
  )
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92
  ))
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$CSR_binary_92
)

#Run ASRs
ASR_AMOnly_Rest_CSR_binary_92_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_92_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_92_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92,
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
    ASR_AMOnly_Rest_CSR_binary_92_ER_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_92_ARD_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_92_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AMOnly_Rest_CSR_binary_92_ER_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_92_ER_yang")
save(ASR_AMOnly_Rest_CSR_binary_92_ARD_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_92_ARD_yang")
save(ASR_AMOnly_Rest_CSR_binary_92_SYM_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_92_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AMOnly_Rest_CSR_binary_92_ARD_yang
plotMKmodel(ASR_AMOnly_Rest_CSR_binary_92_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AMOnly_Rest_CSR_binary_92 <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92 %>%
  dplyr::select(Symbiotic_type_AMOnly_Rest, CSR_binary_92)
row.names(dat_plot_AMOnly_Rest_CSR_binary_92) <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_92$Species_name
head(dat_plot_AMOnly_Rest_CSR_binary_92)


#CSR ASR - plot pdf
pdf("./Output/ASR_AMOnly_Rest_CSR_binary_92.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AMOnly_Rest_CSR_binary_92,
  cols = list(
    Symbiotic_type_AMOnly_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_92 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AMOnly_Rest_CSR_binary_92_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AMOnly_Rest_CSR_binary_92,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AMOnly (AM)","Rest (All other types)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_92$Symbiotic_type_AMOnly_Rest
vec_selection_type_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_92$CSR_binary_92
names(vec_symbiont_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_92)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_92)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AMOnly_Rest_CSR_binary_92_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AMOnly_Rest_CSR_binary_92_ARD
plot(pagel_AMOnly_Rest_CSR_binary_92_ARD)

pdf("./Output/Plot_pagel_AMOnly_Rest_CSR_binary_92_ARD.pdf")
plot(pagel_AMOnly_Rest_CSR_binary_92_ARD)
dev.off()

save.image()

######95

#First we'll study the baseline case (cutoff cosine similarity = 95, AnyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AnyAM_Rest, CSR_binary_95)
head(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95
)
states_AnyAM_Rest_CSR_binary_95 <-
  c("AnyAM & CSR95", "AnyAM & NoCSR", "noAM & CSR95", "noAM & NoCSR")

analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type_AnyAM_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type_AnyAM_Rest
    )
  )
analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95
  ))
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type_AnyAM_Rest)
table(analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95)
table(
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Symbiotic_type,
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$CSR_binary_95
)

#Run ASRs
ASR_AnyAM_Rest_CSR_binary_95_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_95_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AnyAM_Rest_CSR_binary_95_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95,
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
    ASR_AnyAM_Rest_CSR_binary_95_ER_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_95_ARD_yang$AICc,
    ASR_AnyAM_Rest_CSR_binary_95_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AnyAM_Rest_CSR_binary_95_ER_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_95_ER_yang")
save(ASR_AnyAM_Rest_CSR_binary_95_ARD_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_95_ARD_yang")
save(ASR_AnyAM_Rest_CSR_binary_95_SYM_yang, file = "./Output/ASR_AnyAM_Rest_CSR_binary_95_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AnyAM_Rest_CSR_binary_95_ARD_yang
plotMKmodel(ASR_AnyAM_Rest_CSR_binary_95_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AnyAM_Rest_CSR_binary_95 <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95 %>%
  dplyr::select(Symbiotic_type_AnyAM_Rest, CSR_binary_95)
row.names(dat_plot_AnyAM_Rest_CSR_binary_95) <-
  analysis_dat_CSR_symb_AnyAM_Rest_CSR_binary_95$Species_name
head(dat_plot_AnyAM_Rest_CSR_binary_95)

plotvec_symbiont_selection_type_binary_70<-
  c("#e31a1c","#fb9a99",     #Red is AMF
    "#8c510a","#d8b365",     #Brown is AMNod
    "#1f78b4","#a6cee3",      #Blue is ECM
    "#a6cee3","#cab2d6",     #Purple is ECMAM
    "#33a02c","#b2df8a",      #Green is ERM  
    "#fec44f","#fee391",     #Yellow is NM
    "#525252","#bdbdbd",    #Grey is NMAM
    "#014636" ,"#02818a")    #Turqouise is OM 

#CSR ASR - plot pdf
pdf("./Output/ASR_AnyAM_Rest_CSR_binary_95.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AnyAM_Rest_CSR_binary_95,
  cols = list(
    Symbiotic_type_AnyAM_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_95 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AnyAM_Rest_CSR_binary_95_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AnyAM_Rest_CSR_binary_95,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)","Rest (EcM,ErM,NM,OM)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_95$Symbiotic_type_AnyAM_Rest
vec_selection_type_binary <-
  dat_plot_AnyAM_Rest_CSR_binary_95$CSR_binary_95
names(vec_symbiont_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_95)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AnyAM_Rest_CSR_binary_95)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AnyAM_Rest_CSR_binary_95_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AnyAM_Rest_CSR_binary_95_ARD
plot(pagel_AnyAM_Rest_CSR_binary_95_ARD)

pdf("./Output/Plot_pagel_AnyAM_Rest_CSR_binary_95_ARD.pdf")
plot(pagel_AnyAM_Rest_CSR_binary_95_ARD)
dev.off()

save.image()

#Second, the baseline case with OnlyAM (cutoff cosine similarity = 95, OnlyAM vs other)

#Data formatting. We need three columns, species and symbiont state and selection type.
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95 <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type_AMOnly_Rest, CSR_binary_95)
head(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95
)
states_AMOnly_Rest_CSR_binary_95 <-
  c("AMOnly & CSR95", "AMOnly & NoCSR", "noAM & CSR95", "noAM & NoCSR")

analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type_AMOnly_Rest <-
  as.numeric(
    as.factor(
      analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type_AMOnly_Rest
    )
  )
analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95 <-
  as.numeric(as.factor(
    analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95
  ))
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type_AMOnly_Rest)
table(analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95)
table(
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Symbiotic_type,
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$CSR_binary_95
)

#Run ASRs
ASR_AMOnly_Rest_CSR_binary_95_ER_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95,
    rate.cat = 1,
    model = "ER",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_95_ARD_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95,
    rate.cat = 1,
    model = "ARD",
    node.states = "marginal",
    root.p = "yang",
    nstarts = 10,
    n.cores = 7
  )
ASR_AMOnly_Rest_CSR_binary_95_SYM_yang <-
  corHMM(
    phy = analysis_tree,
    data = analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95,
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
    ASR_AMOnly_Rest_CSR_binary_95_ER_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_95_ARD_yang$AICc,
    ASR_AMOnly_Rest_CSR_binary_95_SYM_yang$AICc
  )
)

#Save all model ran
save(ASR_AMOnly_Rest_CSR_binary_95_ER_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_95_ER_yang")
save(ASR_AMOnly_Rest_CSR_binary_95_ARD_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_95_ARD_yang")
save(ASR_AMOnly_Rest_CSR_binary_95_SYM_yang, file = "./Output/ASR_AMOnly_Rest_CSR_binary_95_SYM_yang")

#Let's look at the best (=ARD) model
ASR_AMOnly_Rest_CSR_binary_95_ARD_yang
plotMKmodel(ASR_AMOnly_Rest_CSR_binary_95_ARD_yang)

##Let's for now plot this reconstruction onto the tree

#Create a data frame to plot the trait data
dat_plot_AMOnly_Rest_CSR_binary_95 <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95 %>%
  dplyr::select(Symbiotic_type_AMOnly_Rest, CSR_binary_95)
row.names(dat_plot_AMOnly_Rest_CSR_binary_95) <-
  analysis_dat_CSR_symb_AMOnly_Rest_CSR_binary_95$Species_name
head(dat_plot_AMOnly_Rest_CSR_binary_95)

plotvec_symbiont_selection_type_binary_70<-
  c("#e31a1c","#fb9a99",     #Red is AMF
    "#8c510a","#d8b365",     #Brown is AMNod
    "#1f78b4","#a6cee3",      #Blue is ECM
    "#a6cee3","#cab2d6",     #Purple is ECMAM
    "#33a02c","#b2df8a",      #Green is ERM  
    "#fec44f","#fee391",     #Yellow is NM
    "#525252","#bdbdbd",    #Grey is NMAM
    "#014636" ,"#02818a")    #Turqouise is OM 

#CSR ASR - plot pdf
pdf("./Output/ASR_AMOnly_Rest_CSR_binary_95.pdf",
    width = 20,
    height = 20)
trait.plot(
  tree = analysis_tree,
  dat = dat_plot_AMOnly_Rest_CSR_binary_95,
  cols = list(
    Symbiotic_type_AMOnly_Rest = brewer.pal(n = 8, "Set2"),
    CSR_binary_95 = brewer.pal(n = 3, "Accent")
  ),
  type = "f",
  legend = F,
  w = 1 / 40,
  edge.width = 2,
  cex.lab = 0.01,
  tip.color = "white",
  show.node.label = T
)
nodelabels(pie = ASR_AMOnly_Rest_CSR_binary_95_ARD_yang$states,
           piecol = plotvec_symbiont_binary_selection_type_binary,
           cex = 0.3)
legend(
  legend = states_AMOnly_Rest_CSR_binary_95,
  x = "bottomright",
  fill = plotvec_symbiont_binary_selection_type_binary,
  cex = 1.5
)
legend(
  legend = c("AMOnly (AM)","Rest (All other types)"),
  x = "topleft",
  fill = brewer.pal(n = 8, "Set2"),
  title = "Inner Ring",
  cex=1.5
)
legend(
  legend = c("CSR","NoCSR"),
  x = "topright",
  fill = brewer.pal(n = 3, "Accent"),
  title = "Outer Ring",
  cex=1.5
)
add.scale.bar()
dev.off()

####Run Pagel's model
vec_symbiont_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_95$Symbiotic_type_AMOnly_Rest
vec_selection_type_binary <-
  dat_plot_AMOnly_Rest_CSR_binary_95$CSR_binary_95
names(vec_symbiont_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_95)
names(vec_selection_type_binary) <-
  row.names(dat_plot_AMOnly_Rest_CSR_binary_95)
head(vec_symbiont_binary)
head(vec_selection_type_binary)
table(vec_symbiont_binary)
table(vec_selection_type_binary)

pagel_AMOnly_Rest_CSR_binary_95_ARD <-
  fitPagel(
    tree = analysis_tree,
    x = vec_symbiont_binary,
    y = vec_selection_type_binary,
    model = "ARD",
    pi = "fitzjohn"
  )
pagel_AMOnly_Rest_CSR_binary_95_ARD
plot(pagel_AMOnly_Rest_CSR_binary_95_ARD)

pdf("./Output/Plot_pagel_AMOnly_Rest_CSR_binary_95_ARD.pdf")
plot(pagel_AMOnly_Rest_CSR_binary_95_ARD)
dev.off()

save.image()


