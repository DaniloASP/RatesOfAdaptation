################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#######################################################################################################################
#   This script will summarize GRAPES runs for bootstrap and permutations on the core and accessory categories         #
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
mydata <- read.table("./bootstrap/20spp_Bootstrap.txt", header = TRUE) # this file contains the ouput of 100 bootstrap per species. It was concatenated using cat
mydata <- mydata[!grepl("Species", mydata$Species),] # clean dataframe


#### Data preparation considering the average value
# this file was generated using grep 'Average' *_GrapesResults_allmodels_4categories.txt > 20spp_Average_Grapes.txt
mydata_average <- read.table("./MeanOverModels_grapes/20spp_Average_Grapes.txt") # averaged over all models
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

# plot
#mydata$Species <- factor(mydata$Species, levels=species_vector)
#mydata$Data_f = factor(mydata$Data, levels=c("NonSecreted", "SecretedNonEffector","Effectors"))
#mydata$Data <- gsub("NonSecreted","Non-secreted",mydata$Data)
#mydata$Data <- gsub("SecretedNonEffector","Secreted",mydata$Data)

mydata$Data <- factor(mydata$Data, levels=c("Core", "Accessory"))


#mydata_average$Data <- gsub("Nonsecreted","Non-secreted",mydata_average$Data)
#mydata_average$Data <- gsub("\\<Secreted\\>","SecretedNonEffector",mydata_average$Data)
#mydata_average$Data <- gsub("Effectors","Effectors",mydata_average$Data)
#mydata_average <- mydata_average[!grepl("AllGenes", mydata_average$Data),]
#mydata_average$Data_f = factor(mydata_average$Data, levels=c("NonSecreted", "SecretedNonEffector","Effectors"))



###### common part ends


######################################
# ////// >>>>>> omegaA <<<<<< \\\\\\ #
######################################

####################################### Venturia - VI
# 1
species_vector
Species_name="Venturia"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors



get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_VI <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("V. inaequalis")
p.omegaA_VI

####################################### end for this plot

####################################### Parastagonospora - PN
# 2

Species_name="Parastagonospora"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))



# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_PN <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("P. nodorum")
p.omegaA_PN


####################################### end for this plot

####################################### Sphaerulina - PN
# 3

Species_name="Sphaerulina"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))



# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_SM <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("S. musiva")
p.omegaA_SM


####################################### end for this plot

####################################### Penicillium - PN
# 4 - no data
####################################### end for this plot

####################################### Botrytis - BC
# 5

Species_name="Botrytis"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_BC <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("B. cinera")
p.omegaA_BC


####################################### end for this plot

####################################### Ophiostoma_novo - ON
# 6 - no data
####################################### end for this plot


####################################### Ophiostoma_ulmi - OU
# 7 - no data
####################################### end for this plot

####################################### Zymoseptoria - ZT
# 8
species_vector
Species_name="Zymoseptoria"
Best_model="GammaGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_ZT <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("Z. tritici")
p.omegaA_ZT


####################################### end for this plot

####################################### Ophiostoma_americana - OA
# 9 - no data
####################################### end for this plot

####################################### Magnaporthe_triticum - MOt
# 10 - no data
####################################### end for this plot

####################################### Cercospora - CB
# 11
species_vector
Species_name="Cercospora"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_CB <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("C. beticola")
p.omegaA_CB

####################################### end for this plot

####################################### Aspergilus - AF
# 12
species_vector
Species_name="Aspergilus"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))

# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_AF <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("A. flavus")
p.omegaA_AF

####################################### end for this plot


####################################### Neurospora - ND
# 13
species_vector
Species_name="Neurospora"
Best_model="ScaledBeta"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_ND <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("N. discreta")
p.omegaA_ND

####################################### end for this plot


####################################### Verticillium - VD
# 14 - no data for ACC
####################################### end for this plot

####################################### Sclerotinia - SS
# 15
species_vector
Species_name="Sclerotinia"
Best_model="GammaExpo"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))

# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_SS <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("S. sclerotiorum")
p.omegaA_SS

####################################### end for this plot


####################################### Rhynchosporium - RC
# 16
species_vector
Species_name="Rhynchosporium"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))

# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_RC <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("R. commune")
p.omegaA_RC

####################################### end for this plot

####################################### Pyrenophora - PTT
# 17
species_vector
Species_name="Pyrenophora"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))

# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_PTT <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("P. teres f. teres")
p.omegaA_PTT

####################################### end for this plot



####################################### Magnaporthe_rice - MOr
# 18
species_vector
Species_name="Magnaporthe_rice"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_MOr <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("M. oryzae rice")
p.omegaA_MOr

####################################### end for this plot


####################################### Fusarium - FG
# 19
species_vector
Species_name="Fusarium"
Best_model="DisplGamma"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_FG <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("F. graminearum")
p.omegaA_FG

####################################### end for this plot

####################################### Cryphonectria - CP
# 20
species_vector
Species_name="Cryphonectria"
Best_model="GammaZero"


### Permutation part
path_to_first_run=paste0("../../",{Species_name},"/3_output/3_singleRun")

res.Core.spp <- read.res(paste0({path_to_first_run},"/CORE.genes.res"), "Core")
res.Core.spp[order(res.Core.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Core.spp$Data <- gsub("NonSecreted","Non-secreted",res.Core.spp$Data)
#res.Core.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Core.spp$Data)

res.Accessory.spp <- read.res(paste0({path_to_first_run},"/ACC.genes.res"), "Accessory")
res.Accessory.spp[order(res.Accessory.spp$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]

#res.Accessory.spp$Data <- gsub("NonSecreted","Non-secreted",res.Accessory.spp$Data)
#res.Accessory.spp$Data <- gsub("SecretedNonEffector","Secreted",res.Accessory.spp$Data)

#res.eff.spp <- read.res(paste0({path_to_first_run},"/Effector.genes.res"), "Effectors")
#res.eff.spp[order(res.eff.spp$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA")]#

#res.eff.spp$Data <- gsub("NonSecreted","Non-secreted",res.eff.spp$Data)
#res.eff.spp$Data <- gsub("SecretedNonEffector","Secreted",res.eff.spp$Data)


permutation.all <- read.table(paste0("./permutation/",{Species_name},"_2CoreAcc_Permutation.csv"), header = TRUE, sep = ",")

#permutation.all$Data <- gsub("NonSecreted","Non-secreted",permutation.all$Data)
#permutation.all$Data <- gsub("SecretedNonEffector","Secreted",permutation.all$Data)

perm.Core <- setDT(permutation.all[grepl("Core", permutation.all$Data),]) # perm.clu = perm.Core
perm.Accessory <- setDT(permutation.all[grepl("Accessory", permutation.all$Data),]) # perm.enc = perm.Accessory
#perm.Effectors <- setDT(permutation.all[grepl("Effectors", permutation.all$Data),]) # perm.nec = perm.Effectors

get.perm.pairs <- function(var) {
  perm.pairs <- data.frame(CoreVsAccessory = unlist(perm.Core[,..var] - perm.Accessory[,..var]))
  perm.pairs <- melt(perm.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(perm.pairs)
}

get.dat.pairs <- function(var) {
  dat.pairs <- data.frame(CoreVsAccessory = subset(res.Core.spp, model == {Best_model}, var)[1,1] - subset(res.Accessory.spp, model == {Best_model}, var)[1,1])
  dat.pairs <- melt(dat.pairs, value.name = var, variable.name = "Comparison", measure.vars = 1:1)
  return(dat.pairs)
}

perm.pairs <- get.perm.pairs("omegaA")
dat.pairs <- get.dat.pairs("omegaA")

# Plots with significance levels:

# Test function for getting p-values:
perm.test <- function(var, comp) {
  perm.pairs <- subset(get.perm.pairs(var), Comparison == comp)
  dat.pairs <- subset(get.dat.pairs(var), Comparison == comp)
  return(list(p.value = (sum(abs(unlist(perm.pairs[,var])) >= abs(unlist(dat.pairs[,var])), na.rm = TRUE) + 1) / (sum(!is.na(perm.pairs[,var])) + 1)))
}

format.p.val <- Vectorize(function(p) {
  s <- "(NS)"
  if (p < 0.1) s <- "(.)"
  if (p < 0.05) s <- "(*)"
  if (p < 0.01) s <- "(**)"
  if (p < 0.001) s <- "(***)"
  return(paste(sprintf("%.4f", p), s))
})

format.p.val(p.adjust(perm.test("omegaA", "CoreVsAccessory")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "SecretedVsEffectors")$p.value))
#format.p.val(p.adjust(perm.test("omegaA", "NonSecretedVsEffectors")$p.value))


# plot starts

mydata_species <- mydata[grepl({Species_name}, mydata$Species),]

mydataAverage_test <- mydata_average[grepl({Best_model}, mydata_average$model),]
mydataAverage_test <- mydataAverage_test[grepl({Species_name}, mydataAverage_test$Species),]

p.omegaA_CP <- ggplot(mydata_species, aes(x = Data, y = omegaA)) +
  geom_boxplot() +
  geom_point(data = mydataAverage_test, aes(x = Data, y = omegaA), col = "orange", shape = 18, size = 5) +
  geom_signif(comparisons = list(a = c("Core", "Accessory")),
              annotations = paste("p = ", (format.p.val(c(
                perm.test(var = "omegaA", comp = "CoreVsAccessory")$p.value))), sep = ""), y_position = c((max(mydata_species$omegaA, na.rm = TRUE)+(max(mydata_species$omegaA, na.rm = TRUE)*0.01))), textsize = 2.5, map_signif_level = TRUE) +
  xlab("Category") +
  ylab(expression(omega[A])) +
  theme_bw() + theme(plot.title = element_text(size=10, face="bold.italic"),axis.text.x = element_text(size = 7.5)) + ggtitle("C. parasitica")
p.omegaA_CP

####################################### end for this plot


###########################
# Create panel and export #
###########################

plot_grid(p.omegaA_VI, p.omegaA_PN,p.omegaA_SM, p.omegaA_PB,p.omegaA_BC, p.omegaA_ON,p.omegaA_OU, p.omegaA_ZT,p.omegaA_OA, p.omegaA_MOt)

plot_grid(p.omegaA_VI, p.omegaA_PN,p.omegaA_SM, p.omegaA_BC, p.omegaA_ZT, p.omegaA_CB,p.omegaA_AF,p.omegaA_ND,p.omegaA_SS,p.omegaA_RC,p.omegaA_PTT,p.omegaA_MOr,p.omegaA_FG,p.omegaA_CP, labels = "AUTO", ncol = 4)

ggsave(file = paste0("14spp_omegaA.pdf"), width = 400, height = 400, units = "mm", device='pdf', useDingbats=F)
ggsave(file = paste0("14spp_omegaA.png"), width = 400, height = 400, units = "mm")


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
