################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#######################################################################################################################
#   This script will summarize a GRAPES run on the gene categories using all models available                          #
#                                                                                                                     #
# Files needed:                                                                                                       #
#             * Output files from GRAPES (*.res)                                                                      #
#                                                                                                                     #
#######################################################################################################################
#
#
#
######################################################################################################################
# Created on 14/02/19 by jdutheil
# Plotting modified by Danilo on 12/12/2022


Species_is="Zymoseptoria"
Abb_is="Zymoseptoria"
Abb_slurm="ZT"
Num_chr=485
Best_grapes_model_allgenes="GammaGamma"

# fixed
# libraries and functions
library(data.table)
library(reshape2)
library(plyr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(tidyverse)

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


# ds info for the paper
res.allGenes <- read.res("../../3_output/3_singleRun/AllCoreGenes.genes.res", "AllGenes")
summary_stat_grapes_DF <- res.allGenes[order(res.allGenes$AIC),c("dataset","model", "pn", "Lpn", "ps", "Lps", "dn", "Ldn", "ds", "Lds", "dN","dS","dNdS", "piS", "piN", "pNpS", "omegaNA", "omegaA")]
summary_stat_grapes_DF <- summary_stat_grapes_DF[summary_stat_grapes_DF$model %like% Best_grapes_model_allgenes, ]
summary_stat_grapes_DF$galtier_dS <- (summary_stat_grapes_DF$ds / summary_stat_grapes_DF$dn)
summary_stat_grapes_DF$galtier_dN <- (summary_stat_grapes_DF$Lds / summary_stat_grapes_DF$Ldn)
summary_stat_grapes_DF$galtier_dN_dS <- (summary_stat_grapes_DF$galtier_dN / summary_stat_grapes_DF$galtier_dS)
summary_stat_grapes_DF$species <- Species_is

# export 
summary_stat_grapes_export <- summary_stat_grapes_DF %>% select("species", "pn", "Lpn", "ps", "Lps", "Lds", "Ldn", "ds", "dn", "piS", "piN", "pNpS", "galtier_dS", "galtier_dN", "galtier_dN_dS", "omegaA", "omegaNA")
write.table(summary_stat_grapes_export, paste0(Abb_slurm,"_GrapesResults_PopStat_grapes.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)


#############################################################################################################################################################
### This part was added on the 12 of May of 2023 to extract info on pNpS per gene category for correlating with omegaNA after Julien's feedback on the manuscript ###
#############################################################################################################################################################
# get runs for model selection
res.allGenes <- read.res("../../3_output/3_singleRun/AllCoreGenes.genes.res", "AllGenes")
res.allGenes[order(res.allGenes$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")]

res.NonSecreted <- read.res("../../3_output/3_singleRun/NonSecreted.genes.res", "Nonsecreted")
res.NonSecreted[order(res.NonSecreted$AIC),c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")]

res.Secreted <- read.res("../../3_output/3_singleRun/Secreted.genes.res", "Secreted")
res.Secreted[order(res.Secreted$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")]

res.Effectors <- read.res("../../3_output/3_singleRun/Effector.genes.res", "Effectors")
res.Effectors[order(res.Effectors$AIC), c("model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")]

tbl.FourCategories <- rbind(res.allGenes[, c("Data", "model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")],
                            res.NonSecreted[, c("Data", "model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")],
                            res.Secreted[, c("Data", "model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")],
                            res.Effectors[, c("Data", "model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "pNpS")])

tbl.FourCategories <- subset(tbl.FourCategories, model %in% c("Neutral", "GammaZero", "GammaGamma", "GammaExpo", "DisplGamma", "ScaledBeta"))
# Note: the FGM has an implementation issue, we discard it here.

# Model averaging.
# Computing AIC weights:
m.avg <- function(x, aic) {
  a <- aic
  b <- a - min(a)
  lik <- exp(-0.5*b)
  waic <- lik / sum(lik)
  return(weighted.mean(x, w = waic, na.rm = TRUE))
}

library(plyr)

tbl.FourCategories.avg <- ddply(subset(tbl.FourCategories, model %in% c("Neutral", "GammaZero", "GammaGamma", "GammaExpo", "DisplGamma", "ScaledBeta")), .variables = "Data", .fun = summarize,
                                alpha = m.avg(alpha, AIC),
                                omegaA = m.avg(omegaA, AIC),
                                omegaNA = m.avg(omegaNA, AIC),
                                model = "Average")

tbl.FourCategories.2 <- merge(tbl.FourCategories, tbl.FourCategories.avg, all=TRUE)
tbl.FourCategories.3 <- ddply(tbl.FourCategories.2, .variables = "Data", .fun = function(d) {
  rownames(d) <- d$model
  d[c("Neutral", "GammaZero", "GammaGamma", "GammaExpo", "DisplGamma", "ScaledBeta", "Average"),
    c("Data", "model", "lnL", "AIC", "alpha", "omegaA", "omegaNA")]
})
tbl.FourCategories.export <- tbl.FourCategories.3
tbl.FourCategories.export
tbl.FourCategories.export$lnL <- format(round(tbl.FourCategories.export$lnL, 3), nsmall = 3)
tbl.FourCategories.export$AIC <- format(round(tbl.FourCategories.export$AIC, 3), nsmall = 3)
tbl.FourCategories.export$alpha <- format(round(tbl.FourCategories.export$alpha, 3), nsmall = 3)
tbl.FourCategories.export$omegaA <- format(round(tbl.FourCategories.export$omegaA, digits = 3), nsmall = 3)
tbl.FourCategories.export$omegaNA <- format(round(tbl.FourCategories.export$omegaNA, 3), nsmall = 3)
tbl.FourCategories.export$Species <- Species_is

# add pNpS to final table
tbl.FourCategories.export$merge_ID <- paste0(tbl.FourCategories.export$Data,"_",tbl.FourCategories.export$model)
tbl.FourCategories$merge_ID <- paste0(tbl.FourCategories$Data,"_",tbl.FourCategories$model)

tbl.FourCategories.export <- merge(x = tbl.FourCategories.export, y = tbl.FourCategories[ , c("merge_ID", "pNpS")], by = "merge_ID", all.x=TRUE)

tbl.FourCategories.export$merge_ID <- NULL

tbl.FourCategories.export$pNpS <- format(round(tbl.FourCategories.export$pNpS, 4), nsmall = 4)

# export 
write.table(tbl.FourCategories.export, paste0(Abb_slurm,"_GrapesResults_PopStat_grapes_pNpS-omegaNA.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)



######################################################################################################################
#                                                      END                                                           #
######################################################################################################################


