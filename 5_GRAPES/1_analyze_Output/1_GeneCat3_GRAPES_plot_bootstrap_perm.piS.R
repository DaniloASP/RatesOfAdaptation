################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#######################################################################################################################
#   This script will summarize GRAPES runs for bootstrap and permutations on the gene categories                       #
#                                                                                                                     #
# Files needed:                                                                                                       #
#             * Output files from GRAPES (*.res)                                                                      #
#                                                                                                                     #
#######################################################################################################################
#
#
#
######################################################################################################################


library(data.table)
library(reshape2)
library(plyr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(tidyverse)
library(doBy)
library(ggpubr)
library(patchwork)


###### common part starts


# function
read.res <- function(file, name) {
  res <- read.csv(file, stringsAsFactors = FALSE)
  res$dN <- (res$dn / res$Ldn) 
  res$dS <- (res$ds / res$Lds) 
  res$dNdS <- (res$dn / res$Ldn) / (res$ds / res$Lds) 
  res$dNdS_obs <- (res$obs_dN / res$Ldn) / (res$obs_dS / res$Lds) 
  res$dNdS_exp <- (res$exp_dN / res$Ldn) / (res$exp_dS / res$Lds) 
  res$piN <- (res$pn / res$Lpn)
  res$piS <- (res$ps / res$Lps) 
  res$pNpS <- (res$pn / res$Lpn) / (res$ps / res$Lps) 
  res$omegaNA <- with(res, dNdS - omegaA)
  res$Data <- name
  res$AIC <- 2*res$X.param - 2*res$lnL
  return(res)
}


# read data with bootstrap information
mydata <- read.table("./results/20spp_Bootstrap.txt", header = TRUE) # this file contains the ouput of 100 bootstrap per species. It was concatenated using bash
mydata <- mydata[!grepl("Species", mydata$Species),] # clean dataframe

#### Data preparation considering the average value
# this file was generated using grep 'Average' *_GrapesResults_allmodels_4categories.txt > 20spp_Average_Grapes.txt
mydata_average <- read.table("../MeanOverModels_grapes/20spp_Average_Grapes.txt") # averaged over all models
colnames(mydata_average) <- c("Data","model","lnL","AIC","alpha","omegaA","omegaNA","Species")
mydata_average$Data <- gsub(".*:","",mydata_average$Data)
mydata_average$Data

# the order is based on the average across all models for omegaA *allgenes*
#mydata_average_allgenes <- mydata_average[grepl("AllGenes", mydata_average$Data),]
#mydata_average_allgenes <- mydata_average_allgenes[grepl("Average", mydata_average_allgenes$model),]

#mydata_average_allgenes$Species = with(mydata_average_allgenes, reorder(Species, omegaA))
species_vector <- as.vector(unique(mydata$Species))

### Preparation for plotting
# remove headers in rows
mydata <- mydata[!grepl("Data", mydata$Data),]

# convert str to numeric
mydata$omegaA <- as.numeric(mydata$omegaA)
mydata$omegaNA <- as.numeric(mydata$omegaNA)
mydata$alpha <- as.numeric(mydata$alpha)
mydata$piN <- as.numeric(mydata$piN)
mydata$piS <- as.numeric(mydata$piS)
mydata$pNpS <- as.numeric(mydata$pNpS)
mydata$dN <- as.numeric(mydata$dN)
mydata$dS <- as.numeric(mydata$dS)
mydata$dNdS <- as.numeric(mydata$dNdS)

# plot
#mydata$Species <- factor(mydata$Species, levels=species_vector)
#mydata$Data_f = factor(mydata$Data, levels=c("NonSecreted", "SecretedNonEffector","Effectors"))
mydata$Data <- gsub("NonSecreted","Non-secreted",mydata$Data)
mydata$Data <- gsub("SecretedNonEffector","Secreted",mydata$Data)

mydata$Data <- factor(mydata$Data, levels=c("Non-secreted", "Secreted", "Effectors"))


#mydata_average$Data <- gsub("Nonsecreted","Non-secreted",mydata_average$Data)
#mydata_average$Data <- gsub("\\<Secreted\\>","SecretedNonEffector",mydata_average$Data)
#mydata_average$Data <- gsub("Effectors","Effectors",mydata_average$Data)
#mydata_average <- mydata_average[!grepl("AllGenes", mydata_average$Data),]
#mydata_average$Data_f = factor(mydata_average$Data, levels=c("NonSecreted", "SecretedNonEffector","Effectors"))



###### common part ends

# Note that the variable here can be changed to any of the other calculations, e.g, omega, omegaNA etc
######################################
# ////// >>>>>> piS <<<<<< \\\\\\ #
######################################

####################################### Venturia - VI
# 1
species_vector
Species_name="Venturia"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors



get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsSecreted")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_VI <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS)+(max(mydata_species$piS)*0.01)), (max(mydata_species$piS)+(max(mydata_species$piS)*0.1)), (max(mydata_species$piS)+(max(mydata_species$piS)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("V. inaequalis")
p.piS_VI

####################################### end for this plot

####################################### Parastagonospora - PN
# 2

Species_name="Parastagonospora"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_PN <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("P. nodorum")
p.piS_PN


####################################### end for this plot

####################################### Sphaerulina - PN
# 3

Species_name="Sphaerulina"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_SM <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("S. musiva")
p.piS_SM


####################################### end for this plot

####################################### Penicillium - PN
# 4

Species_name="Penicillium"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_PB <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("P. biforme")
p.piS_PB


####################################### end for this plot

####################################### Botrytis - BC
# 5

Species_name="Botrytis"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_BC <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("B. cinera")
p.piS_BC


####################################### end for this plot

####################################### Ophiostoma_novo - ON
# 6
species_vector
Species_name="Ophiostoma_novo"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_ON <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("O. novo")
p.piS_ON


####################################### end for this plot


####################################### Ophiostoma_ulmi - OU
# 7
species_vector
Species_name="Ophiostoma_ulmi"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_OU <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("O. ulmi")
p.piS_OU


####################################### end for this plot

####################################### Zymoseptoria - ZT
# 8
species_vector
Species_name="Zymoseptoria"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_ZT <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("Z. tritici")
p.piS_ZT


####################################### end for this plot

####################################### Ophiostoma_americana - OA
# 9
species_vector
Species_name="Ophiostoma_americana"
Best_model="ScaledBeta"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_OA <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("O. americana")
p.piS_OA


####################################### end for this plot

####################################### Magnaporthe_triticum - MOt
# 10
species_vector
Species_name="Magnaporthe_triticum"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_MOt <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("M. oryzae Triticum")
p.piS_MOt


####################################### end for this plot

####################################### Cercospora - CB
# 11
species_vector
Species_name="Cercospora"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_CB <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("C. beticola")
p.piS_CB

####################################### end for this plot

####################################### Aspergilus - AF
# 12
species_vector
Species_name="Aspergilus"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_AF <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("A. flavus")
p.piS_AF

####################################### end for this plot


####################################### Neurospora - ND
# 13
species_vector
Species_name="Neurospora"
Best_model="ScaledBeta"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_ND <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("N. discreta")
p.piS_ND

####################################### end for this plot


####################################### Verticillium - VD
# 14
species_vector
Species_name="Verticillium"
Best_model="GammaZero"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_VD <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("V. dahliae")
p.piS_VD

####################################### end for this plot

####################################### Sclerotinia - SS
# 15
species_vector
Species_name="Sclerotinia"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_SS <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("S. sclerotiorum")
p.piS_SS

####################################### end for this plot


####################################### Rhynchosporium - RC
# 16
species_vector
Species_name="Rhynchosporium"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_RC <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("R. commune")
p.piS_RC

####################################### end for this plot

####################################### Pyrenophora - PTT
# 17
species_vector
Species_name="Pyrenophora"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_PTT <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("P. teres f. teres")
p.piS_PTT

####################################### end for this plot



####################################### Magnaporthe_rice - MOr
# 18
species_vector
Species_name="Magnaporthe_rice"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_MOr <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("M. oryzae rice")
p.piS_MOr

####################################### end for this plot


####################################### Fusarium - FG
# 19
species_vector
Species_name="Fusarium"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_FG <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("F. graminearum")
p.piS_FG

####################################### end for this plot

####################################### Cryphonectria - CP
# 20
species_vector
Species_name="Cryphonectria"
Best_model="GammaZero"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.allG.spp <- read.res(paste0({path_to_first_run},"/NonSecreted.genes.res"), "NonSecreted")
res.allG.spp[order(res.allG.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.allG.spp$Data <- gsub("NonSecreted","Non-secreted",res.allG.spp$Data)
res.allG.spp$Data <- gsub("SecretedNonEffector","Secreted",res.allG.spp$Data)

res.NonEff.spp <- read.res(paste0({path_to_first_run},"/Secreted.genes.res"), "SecretedNonEffector")
res.NonEff.spp[order(res.NonEff.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.NonEff.spp$Data <- gsub("NonSecreted","Non-secreted",res.NonEff.spp$Data)
res.NonEff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.NonEff.spp$Data)

res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_3GeneCategories_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.NonSecreted <- setDT(permutation.all[grepl("Non-secreted", permutation.all$Data),]) # perm.clu = perm.NonSecreted
perm.Secreted <- setDT(permutation.all[grepl("Secreted", permutation.all$Data),]) # perm.enc = perm.Secreted
perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(NonSecretedVsSecreted = unlist(perm.NonSecreted[,..var] - perm.Secreted[,..var]),
                           SecretedVsEffectors = unlist(perm.Secreted[,..var] - perm.Effectors[,..var]),
                           NonSecretedVsEffectors = unlist(perm.NonSecreted[,..var] - perm.Effectors[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(NonSecretedVsSecreted = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.NonEff.spp, model == {Best_model}, var)[1,1], SecretedVsEffectors = subset(res.NonEff.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1], NonSecretedVsEffectors = subset(res.allG.spp, model == {Best_model}, var)[1,1] - subset(res.eff.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:3)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("piS")
dat.pairs <- get.dat.pairs("piS")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}
perm.test("omegaA", "NonSecretedVsSecreted")$p.value
perm.test("omegaA", "SecretedVsEffectors")$p.value
perm.test("omegaA", "NonSecretedVsEffectors")$p.value


format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.piS_CP <- ggplot(mydata_species, aes(x = Data, y = piS)) +
  geom_boxplot() +
  geom_signif(comparisons = list(a = c("Non-secreted", "Secreted"),
                                 b = c("Secreted", "Effectors"),
                                 c = c("Non-secreted", "Effectors")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "piS", comp = "NonSecretedVsSecreted")$p.value,
                perm.test(var = "piS", comp = "SecretedVsEffectors")$p.value,
                perm.test(var = "piS", comp = "NonSecretedVsEffectors")$p.value))), sep = ""), y_position = c((max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.01)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.1)), (max(mydata_species$piS, na.rm = TRUE)+(max(mydata_species$piS, na.rm = TRUE)*0.18))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(pi[S])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5),strip.background = element_rect(fill="orange")) + ggtitle("C. parasitica")
p.piS_CP

####################################### end for this plot


###########################
# Create panel and export #
###########################

plot_grid(p.piS_VI, p.piS_PN,p.piS_SM, p.piS_PB,p.piS_BC, p.piS_ON,p.piS_OU, p.piS_ZT,p.piS_OA, p.piS_MOt)

plot_grid(p.piS_VI, p.piS_PN,p.piS_SM, p.piS_PB,p.piS_BC, p.piS_ON,p.piS_OU, p.piS_ZT,p.piS_OA, p.piS_MOt,p.piS_CB,p.piS_AF,p.piS_ND,p.piS_SS,p.piS_VD,p.piS_RC,p.piS_PTT,p.piS_MOr,p.piS_FG,p.piS_CP, labels = "AUTO", ncol = 4)

#ggsave(file = paste0("20spp_piS.pdf"), width = 400, height = 350, units = "mm", device='pdf', useDingbats=F)
ggsave(file = paste0("20spp_piS.pdf"), width = 400, height = 400, units = "mm", device='pdf', useDingbats=F)
ggsave(file = paste0("20spp_piS.png"), width = 400, height = 400, units = "mm")


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################


