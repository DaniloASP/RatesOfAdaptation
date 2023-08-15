################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#######################################################################################################################
#   This script will summarize a GRAPES run on the recombination categories using all models available                #
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
# Plotting modified by Danilo on 30/05/2021


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
  res$dN <- (res$Lds / res$Ldn) ##### modified
  res$dS <- (res$ds / res$dn) ##### modified
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

# get runs for model selection
list_of_files <- system("ls ../../2_output/*.res", intern = TRUE)
#test_file <- list_of_files[1]

Final_DF <- data_frame() # create an empty dataframe
for (output_file in list_of_files) {
  
  recombination_category <- basename(output_file)
  recombination_category <- gsub(".*_D","", recombination_category)
  recombination_category <- gsub(".res","", recombination_category)
  
  res.allGenes <- read.res(output_file, paste0(recombination_category)) # paste the recombination value
  res.allGenes <- res.allGenes[order(res.allGenes$AIC),c("Data","model", "X.param","lnL", "AIC", "alpha", "omegaA", "omegaNA", "piS")]
  res.allGenes <- res.allGenes[grepl(paste0(Best_grapes_model_allgenes), res.allGenes$model),]
  print(res.allGenes)
  
  Final_DF <- rbind(res.allGenes,Final_DF)
}

Final_DF$Species <- Species_is

write.table(Final_DF, paste0(Abb_slurm,"_GrapesResults_RECcategories+piS.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################


