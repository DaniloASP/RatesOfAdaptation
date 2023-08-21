################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
################################################################################################################
#   This script will summarize the output of ANGSD into a single table                                         #
#                                                                                                              #
# Files needed:                                                                                                #
#             * output files for ANGSD calculation of Watterson's estimator theta                              #
#                                                                                                              #
################################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################


# Get output from ANGSD for all species and make a single table
#
#libs
library(dplyr)
library(stringr) # https://stackoverflow.com/questions/39086400/extracting-a-string-between-other-two-strings-in-r

#
output_species_chr <- system("ls *.pestPG", intern = TRUE)

# for testing
#output_species_chr <-output_species_chr[3]

												
# loop
for (species in output_species_chr) {
  #print(species)
  # load data
  wTheta_DF <- read.table(species, header = FALSE)
  
  # rename columns
  colnames(wTheta_DF) <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
  
  # get spp ID
  spp_ID <- str_match(species, "(.*?)\\s*.theta")
  spp_ID <- spp_ID[,2]
  
  # get watterson's theta per bp
  wTheta_DF$tW_bp <- wTheta_DF$tW / wTheta_DF$nSites
  #print(median(wTheta_DF$tW_bp, na.rm = TRUE))
  #print(mean(wTheta_DF$tW_bp, na.rm = TRUE))
  mean_tW_bp <- format(round(mean(wTheta_DF$tW_bp, na.rm = TRUE), 5), nsmall = 5)
  median_tW_bp <- format(round(median(wTheta_DF$tW_bp, na.rm = TRUE), 5), nsmall = 5)
  
  # make the df
  data.WattTheta <- data.frame(spp_ID, mean_tW_bp, median_tW_bp)
  
  ### write full data 
  system("mkdir output_tw")
  write.table(data.WattTheta,paste0("output_tw/",spp_ID,".tW.txt"), quote=F, sep="\t", row.names = F, col.names=F)
  
}


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
