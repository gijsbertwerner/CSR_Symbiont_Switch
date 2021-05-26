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
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_dif_cutoffs_01-12-2020_forR.csv",
           as.is = T,
           strip.white = T)
  # dat_CSR_symb <-
  #   dat_CSR_symb[sample(1:nrow(dat_CSR_symb), size = 100), ] #Only when testing the script
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
#All present

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

#Write the analysed tree
write.tree(phy = analysis_tree,file = "./Output/Cosme_Analysed_Tree.newick")
write.nexus(phy = analysis_tree,file = "./Output/Cosme_Analysed_Tree.nexus")

#### Match dataset and phylogeny

analysis_dat_CSR_symb <-
  dat_CSR_symb %>% filter(Species_name %in% analysis_tree$tip.label)
nrow(analysis_dat_CSR_symb)
#More species in data set, than in tree.


# Analysing the symbiotic states ------------------------------------------

#Data formatting. We need two columns, species and symbiont state.
analysis_dat_CSR_symb_ASR_symbiont_type <-
  analysis_dat_CSR_symb %>% dplyr::select(Species_name, Symbiotic_type)
head(analysis_dat_CSR_symb_ASR_symbiont_type)
table(analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type)
analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type <-
  as.numeric(as.factor(analysis_dat_CSR_symb_ASR_symbiont_type$Symbiotic_type))
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
    nstarts = 25,
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
    nstarts = 25,
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
    nstarts = 25,
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
legend(
  legend = names(table(analysis_dat_CSR_symb$Symbiotic_type)),
  x = "bottomright",
  fill = brewer.pal(n = 8, "Set2"),
  cex = 2
)
add.scale.bar()
dev.off()


# Correlated evolution between the two variables --------------------------

####Colouring for plotting later on. 
plotvec_symbiont_binary_selection_type_binary <-
  c(
    "#e31a1c",
    "#fb9a99",
    #Red is AMF
    "#8c510a",
    "#d8b365",
    #Brown is AMNod
    "#1f78b4",
    "#a6cee3",
    #Blue is ECM
    "#a6cee3",
    "#cab2d6",
    #Purple is ECMAM
    "#33a02c",
    "#b2df8a",
    #Green is ERM
    "#fec44f",
    "#fee391",
    #Yellow is NM
    "#525252",
    "#bdbdbd",
    #Grey is NMAM
    "#014636" ,
    "#02818a"
  )    #Turqouise is OM


### Loading data  for correlated evolution 

# Correlated AnyAMvsnoAM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AnyAMvsnoAM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AnyAMvsnoAMtab.csv",
           as.is = T,
           strip.white = T)
 # dat_CSR_symb_AnyAMvsnoAM <-
 #   dat_CSR_symb_AnyAMvsnoAM[sample(1:nrow(dat_CSR_symb_AnyAMvsnoAM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AnyAMvsnoAM)
table(dat_CSR_symb_AnyAMvsnoAM$Binary_Symb_anyAMvsnoAM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AnyAMvsnoAM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AnyAMvsnoAM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AnyAMvsnoAM$Species_name)
  )
analysis_tree_AnyAMvsnoAM

###General data prep

head(dat_CSR_symb_AnyAMvsnoAM)
analysis_dat_CSR_symb_AnyAMvsnoAM <-
  dat_CSR_symb_AnyAMvsnoAM %>% dplyr::select(-Symbiotic_type,
                                             -C.selection,
                                             -S.selection,
                                             -R.selection,
                                             -Cosine_CSR)
head(analysis_dat_CSR_symb_AnyAMvsnoAM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AnyAMvsnoAM)
states_print_label <- c("AnyAM & CSR_gen",
                        "AnyAM & CSR_spe",
                        "noAM & CSR_gen",
                        "noAM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AnyAMvsnoAM

for (i in 1:(ncol(analysis_dat_CSR_symb_AnyAMvsnoAM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AnyAMvsnoAM[, c(1, 2, i + 2)]
  
  print(paste("This is AnyAMvsnoAM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AnyAMvsnoAM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AnyAMvsnoAM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_anyAMvsnoAM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AnyAM (AM,NMAM,AMNod,EcMAM)", "NoAM (EcM,ErM,NM,OM)"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AnyAMvsnoAM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}

# Correlated AMvsECM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AMvsECM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_onlyAMvsonlyECMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AMvsECM <-
#   dat_CSR_symb_AMvsECM[sample(1:nrow(dat_CSR_symb_AMvsECM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AMvsECM)
table(dat_CSR_symb_AMvsECM$Binary)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AMvsECM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AMvsECM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AMvsECM$Species_name)
  )
analysis_tree_AMvsECM

###General data prep

head(dat_CSR_symb_AMvsECM)
analysis_dat_CSR_symb_AMvsECM <-
  dat_CSR_symb_AMvsECM %>% dplyr::select(    -C.selection,
                                             -S.selection,
                                             -R.selection,
                                             -Cosine_CSR)
head(analysis_dat_CSR_symb_AMvsECM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AMvsECM)
states_print_label <- c("onlyAM & CSR_gen",
                        "onlyAM & CSR_spe",
                        "onlyECM & CSR_gen",
                        "onlyECM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AMvsECM

for (i in 1:(ncol(analysis_dat_CSR_symb_AMvsECM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AMvsECM[, c(1, 2, i + 2)]
  
  print(paste("This is onlyAMvsonlyECM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_onlyAMvsonlyECM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_onlyAMvsonlyECM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AMvsECM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM", "ECM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_onlyAMvsonlyECM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}





# Correlated AM_plus_EcMAMvsEcM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AM_plus_EcMAMvsEcM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AM_plus_EcMAMvsEcMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AM_plus_EcMAMvsEcM <-
 #  dat_CSR_symb_AM_plus_EcMAMvsEcM[sample(1:nrow(dat_CSR_symb_AM_plus_EcMAMvsEcM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AM_plus_EcMAMvsEcM)
table(dat_CSR_symb_AM_plus_EcMAMvsEcM$Binary_Symb_AM_plus_EcMAMvsEcM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AM_plus_EcMAMvsEcM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AM_plus_EcMAMvsEcM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AM_plus_EcMAMvsEcM$Species_name)
  )
analysis_tree_AM_plus_EcMAMvsEcM

###General data prep

head(dat_CSR_symb_AM_plus_EcMAMvsEcM)
analysis_dat_CSR_symb_AM_plus_EcMAMvsEcM <-
  dat_CSR_symb_AM_plus_EcMAMvsEcM %>% dplyr::select(-Symbiotic_type,
                                             -C.selection,
                                             -S.selection,
                                             -R.selection,
                                             -Cosine_CSR)
head(analysis_dat_CSR_symb_AM_plus_EcMAMvsEcM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AM_plus_EcMAMvsEcM)
states_print_label <- c("AM_plus_EcMAM & CSR_gen",
                        "AM_plus_EcMAM & CSR_spe",
                        "EcM & CSR_gen",
                        "EcM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AM_plus_EcMAMvsEcM

for (i in 1:(ncol(analysis_dat_CSR_symb_AM_plus_EcMAMvsEcM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AM_plus_EcMAMvsEcM[, c(1, 2, i + 2)]
  
  print(paste("This is AM_plus_EcMAMvsEcM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AM_plus_EcMAMvsEcM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AM_plus_EcMAMvsEcM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AM_plus_EcMAMvsEcM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM_plus_EcMAM", "EcM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AM_plus_EcMAMvsEcM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}


# Correlated AMvsEcM_plus_EcMAM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AMvsEcM_plus_EcMAM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AMvsEcM_plus_EcMAMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AMvsEcM_plus_EcMAM <-
#  dat_CSR_symb_AMvsEcM_plus_EcMAM[sample(1:nrow(dat_CSR_symb_AMvsEcM_plus_EcMAM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AMvsEcM_plus_EcMAM)
table(dat_CSR_symb_AMvsEcM_plus_EcMAM$Binary_Symb_AMvsEcM_plus_EcMAM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AMvsEcM_plus_EcMAM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AMvsEcM_plus_EcMAM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AMvsEcM_plus_EcMAM$Species_name)
  )
analysis_tree_AMvsEcM_plus_EcMAM

###General data prep

head(dat_CSR_symb_AMvsEcM_plus_EcMAM)
analysis_dat_CSR_symb_AMvsEcM_plus_EcMAM <-
  dat_CSR_symb_AMvsEcM_plus_EcMAM %>% dplyr::select(-Symbiotic_type,
                                                    -C.selection,
                                                    -S.selection,
                                                    -R.selection,
                                                    -Cosine_CSR)
head(analysis_dat_CSR_symb_AMvsEcM_plus_EcMAM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AMvsEcM_plus_EcMAM)
states_print_label <- c("AM & CSR_gen",
                        "AM & CSR_spe",
                        "EcM_plus_EcMAM & CSR_gen",
                        "EcM_plus_EcMAM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AMvsEcM_plus_EcMAM

for (i in 1:(ncol(analysis_dat_CSR_symb_AMvsEcM_plus_EcMAM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AMvsEcM_plus_EcMAM[, c(1, 2, i + 2)]
  
  print(paste("This is AMvsEcM_plus_EcMAM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AMvsEcM_plus_EcMAM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AMvsEcM_plus_EcMAM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AMvsEcM_plus_EcMAM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM", "EcM_plus_EcMAM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AMvsEcM_plus_EcMAM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}


# Correlated AMvsNM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AMvsNM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AMvsNMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AMvsNM <-
#  dat_CSR_symb_AMvsNM[sample(1:nrow(dat_CSR_symb_AMvsNM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AMvsNM)
table(dat_CSR_symb_AMvsNM$Binary_Symb_AMvsNM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AMvsNM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AMvsNM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AMvsNM$Species_name)
  )
analysis_tree_AMvsNM

###General data prep

head(dat_CSR_symb_AMvsNM)
analysis_dat_CSR_symb_AMvsNM <-
  dat_CSR_symb_AMvsNM %>% dplyr::select(            -C.selection,
                                                    -S.selection,
                                                    -R.selection,
                                                    -Cosine_CSR)
head(analysis_dat_CSR_symb_AMvsNM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AMvsNM)
states_print_label <- c("AM & CSR_gen",
                        "AM & CSR_spe",
                        "NM & CSR_gen",
                        "NM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AMvsNM

for (i in 1:(ncol(analysis_dat_CSR_symb_AMvsNM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AMvsNM[, c(1, 2, i + 2)]
  
  print(paste("This is AMvsNM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AMvsNM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AMvsNM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AMvsNM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM", "NM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AMvsNM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}


# Correlated AM_plus_NMAMvsNM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AM_plus_NMAMvsNM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AM_plus_NMAMvsNMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AM_plus_NMAMvsNM <-
#  dat_CSR_symb_AM_plus_NMAMvsNM[sample(1:nrow(dat_CSR_symb_AM_plus_NMAMvsNM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AM_plus_NMAMvsNM)
table(dat_CSR_symb_AM_plus_NMAMvsNM$Binary_Symb_AM_plus_NMAMvsNM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AM_plus_NMAMvsNM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AM_plus_NMAMvsNM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AM_plus_NMAMvsNM$Species_name)
  )
analysis_tree_AM_plus_NMAMvsNM

###General data prep

head(dat_CSR_symb_AM_plus_NMAMvsNM)
analysis_dat_CSR_symb_AM_plus_NMAMvsNM <-
  dat_CSR_symb_AM_plus_NMAMvsNM %>% dplyr::select(-Symbiotic_type,
                                        -C.selection,
                                        -S.selection,
                                        -R.selection,
                                        -Cosine_CSR)
head(analysis_dat_CSR_symb_AM_plus_NMAMvsNM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AM_plus_NMAMvsNM)
states_print_label <- c("AM_plus_NMAM & CSR_gen",
                        "AM_plus_NMAM & CSR_spe",
                        "NM & CSR_gen",
                        "NM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AM_plus_NMAMvsNM

for (i in 1:(ncol(analysis_dat_CSR_symb_AM_plus_NMAMvsNM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AM_plus_NMAMvsNM[, c(1, 2, i + 2)]
  
  print(paste("This is AM_plus_NMAMvsNM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AM_plus_NMAMvsNM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AM_plus_NMAMvsNM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AM_plus_NMAMvsNM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM_plus_NMAM", "NM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AM_plus_NMAMvsNM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}

# Correlated AMvsNM_plus_NMAM -----------------------------------------------------------

#Loading data
dat_CSR_symb_AMvsNM_plus_NMAM <-
  read.csv(file = "./Data/AnalysedData_Gijsbert_23-11-2020_CSR_generalist_NoAM_EcM_NM_comparisons_19-03-2021_AMvsNM_plus_NMAMtab.csv",
           as.is = T,
           strip.white = T)
# dat_CSR_symb_AMvsNM_plus_NMAM <-
#  dat_CSR_symb_AMvsNM_plus_NMAM[sample(1:nrow(dat_CSR_symb_AMvsNM_plus_NMAM), size = 25), ] #Only when testing the script
head(dat_CSR_symb_AMvsNM_plus_NMAM)
table(dat_CSR_symb_AMvsNM_plus_NMAM$Binary_Symb_AMvsNM_plus_NMAM)

####Clean data file
#How many of the species in the database are absent in Smith and Brown?
length(setdiff(dat_CSR_symb_AMvsNM_plus_NMAM$Species_name, smith_brown_tree$tip.label))
#All present

###Clean phylogeny
analysis_tree_AMvsNM_plus_NMAM <-
  drop.tip(
    phy = smith_brown_tree,
    tip = setdiff(smith_brown_tree$tip.label, dat_CSR_symb_AMvsNM_plus_NMAM$Species_name)
  )
analysis_tree_AMvsNM_plus_NMAM

###General data prep

head(dat_CSR_symb_AMvsNM_plus_NMAM)
analysis_dat_CSR_symb_AMvsNM_plus_NMAM <-
  dat_CSR_symb_AMvsNM_plus_NMAM %>% dplyr::select(-Symbiotic_type,
                                                  -C.selection,
                                                  -S.selection,
                                                  -R.selection,
                                                  -Cosine_CSR)
head(analysis_dat_CSR_symb_AMvsNM_plus_NMAM)

#Analysis

###Analysis run specific prepping. 
ncol(analysis_dat_CSR_symb_AMvsNM_plus_NMAM)
states_print_label <- c("AM & CSR_gen",
                        "AM & CSR_spe",
                        "NM_plus_NMAM & CSR_gen",
                        "NM_plus_NMAM & CSR_spe")
analysis_tree_correlate_run<-analysis_tree_AMvsNM_plus_NMAM

for (i in 1:(ncol(analysis_dat_CSR_symb_AMvsNM_plus_NMAM) - 2)) {
  #Data formatting. We need three columns, speies and symbiont state and selection type.
  analysis_dat_frame <-
    analysis_dat_CSR_symb_AMvsNM_plus_NMAM[, c(1, 2, i + 2)]
  
  print(paste("This is AMvsNM_plus_NMAM run:", colnames(analysis_dat_frame)[3]))
  print(paste("It's currently:", Sys.time()))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  analysis_dat_frame[, 2] <-
    as.numeric(as.factor(analysis_dat_frame[, 2]))
  analysis_dat_frame[, 3] <-
    as.numeric(as.factor(analysis_dat_frame[, 3]))
  
  print(table(analysis_dat_frame[, 2],
              analysis_dat_frame[, 3]))
  
  #Run and save the model
  ASR_run_ARD <-
    corHMM(
      phy = analysis_tree_correlate_run,
      data = analysis_dat_frame,
      rate.cat = 1,
      model = "ARD",
      node.states = "marginal",
      root.p = "yang",
      nstarts = 25,
      n.cores = 7
    )
  
  save(ASR_run_ARD,
       file = paste0(
         "./Output/ASR_AMvsNM_plus_NMAM_",
         colnames(analysis_dat_frame)[3],
         "_ARD_yang"
       ))
  
  #Plot the model run
  plotMKmodel(ASR_run_ARD)
  ASR_run_ARD
  
  #Create a data frame to plot the trait data
  dat_plot <-
    analysis_dat_frame[, c(2, 3)]
  row.names(dat_plot) <-
    analysis_dat_frame$Species_name
  colnames(dat_plot)[2] <- "CSR_var"
  
  #CSR ASR - plot pdf
  pdf(
    paste0(
      "./Output/ASR_AMvsNM_plus_NMAM_",
      colnames(analysis_dat_frame)[3],
      ".pdf"
    ),
    width = 20,
    height = 20
  )
  trait.plot(
    tree = analysis_tree_correlate_run,
    dat = dat_plot,
    cols = list(
      Binary_Symb_AMvsNM_plus_NMAM = brewer.pal(n = 8, "Set2"),
      CSR_var = brewer.pal(n = 3, "Accent")
    ),
    type = "f",
    legend = F,
    w = 1 / 40,
    edge.width = 2,
    cex.lab = 0.01,
    tip.color = "white",
    show.node.label = T
  )
  nodelabels(pie = ASR_run_ARD$states,
             piecol = plotvec_symbiont_binary_selection_type_binary,
             cex = 0.3)
  legend(
    legend = states_print_label,
    x = "bottomright",
    fill = plotvec_symbiont_binary_selection_type_binary,
    cex = 1.5
  )
  legend(
    legend = c("AM", "NM_plus_NMAM"),
    x = "topleft",
    fill = brewer.pal(n = 8, "Set2"),
    title = "Inner Ring",
    cex = 1.5
  )
  legend(
    legend = c("CSR_gen", "CSR_spec"),
    x = "topright",
    fill = brewer.pal(n = 3, "Accent"),
    title = "Outer Ring",
    cex = 1.5
  )
  add.scale.bar()
  dev.off()
  
  ####Run Pagel's model
  vec_symbiont_binary <-
    dat_plot[, 1]
  vec_selection_type_binary <-
    dat_plot[, 2]
  names(vec_symbiont_binary) <-
    row.names(dat_plot)
  names(vec_selection_type_binary) <-
    row.names(dat_plot)
  table(vec_symbiont_binary)
  table(vec_selection_type_binary)
  
  pagel_run <-
    fitPagel(
      tree = analysis_tree_correlate_run,
      x = vec_symbiont_binary,
      y = vec_selection_type_binary,
      model = "ARD",
      pi = "fitzjohn"
    )
  pagel_run
  plot(pagel_run)
  
  pdf(paste0(
    "./Output/Plot_pagel_AMvsNM_plus_NMAM",
    colnames(analysis_dat_frame)[3],
    "_ARD.pdf"
  ))
  plot(pagel_run)
  plot.new()
  text(0, 1, paste("likelihood-ratio: ",
                   pagel_run$lik.ratio,
                   "p-value: ",
                   pagel_run$P,
                   collapse='\r\n'), adj = c(0,1), family = 'mono',cex=0.5)
  box()
  dev.off()
  
  #Final cleaning
  gc()
  
}
