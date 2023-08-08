# Created on 13/02/19 by jdutheil
# Modified on 15/05/21 by pereira
# Read all outputs from bpppopstats and reconstruct SFS
# Generate input files for Grapes. Part 1 SFS and output one per gene category, part 2 bootstrap, part 3 permutation

###
# variables according to species

Species_is="Verticillium"
Abb_is="Vdahliae"

Number_of_Isolates=105

# cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/14_grapes/1_first_run
# nano 14_parsePopStatsOutput.R
# module load R/4.1.2
# Rscript 14_parsePopStatsOutput.R
# 
##
input_folder=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/15_popstat/3_fixed_kappa/")

Effector_list_file=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/11_interproscan/effectors.namesONLY.txt")

# modification is only needed to adjust kappa plot
working_Dir_path=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/0_scripts/14_grapes/1_first_run")
setwd(working_Dir_path)

##########
# part 1
##########

###
### Unfolded Site Frequency Spectrum
###

get.sfs.grapes.unfolded <- function(genes, nbChr, prefix) {
  n <- length(genes)
  x.name <- character(n)
  x.nchr <- rep(nbChr, n)
  x.lapo <- numeric(n)
  x.lspo <- numeric(n)
  m.sfsa <- matrix(nrow = n, ncol = nbChr - 1)
  m.sfss <- matrix(nrow = n, ncol = nbChr - 1)
  x.ladi <- numeric(n)
  x.lsdi <- numeric(n)
  x.diva <- numeric(n)
  x.divs <- numeric(n)
  pb <- txtProgressBar(1, n, style = 3)
  for (i in 1:n) {
    setTxtProgressBar(pb, i)
    f <- genes[i]
    x.name[i] <- f
    
    # Get persite calculations:
    sites <- read.table(paste(prefix, f, ".codon.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    sites.com <- subset(sites, MissingDataFrequency == 0) # Only complete sites
    sites.pol <- subset(sites.com, NbAlleles == 2) # Only biallelic sites
    x.lapo[i] <- sum(3 - sites.com$MeanNumberSynPos)
    x.lspo[i] <- sum(sites.com$MeanNumberSynPos)
    sites.pol <- subset(sites.pol, MinorAllele == OutgroupAllele | MajorAllele == OutgroupAllele) # Only sites for which the outgroup allele is one of the ingroup
    sites.pol$AlleleFrequency <- ifelse(sites.pol$MinorAllele == sites.pol$AncestralAllele, sites.pol$MajorAlleleFrequency, sites.pol$MinorAlleleFrequency)
    # Now get the SFSs:
    sitesS <- subset(sites.pol, IsSynPoly == 1)
    sitesA <- subset(sites.pol, IsSynPoly == 0)
    for (j in 1:(nbChr - 1)) {
      m.sfsa[i, j] <- sum(sitesA$AlleleFrequency == j)
      m.sfss[i, j] <- sum(sitesS$AlleleFrequency == j)
    }
    # Divergence patterns:
    #sites.div <- subset(sites, dN + dS <= 1) # Only one mutation per codon
    x.ladi[i] <- sum(3 - sites$MeanNumberSynPosDiv)
    x.lsdi[i] <- sum(sites$MeanNumberSynPosDiv)
    x.diva[i] <- sum(sites$dN)
    x.divs[i] <- sum(sites$dS)
  }
  dat <- cbind(x.lapo, m.sfsa, x.lspo, m.sfss, x.ladi, x.lsdi, x.diva, x.divs)
  return(dat)
}

write.sfs.grapes.unfolded <- function(sfs, file, nbChr, title, dataset) {
  # Write result in Grapes format:
  cat(title, "\n#unfolded\n", dataset, "\t", nbChr, "\t", sep = "", file = file)
  cat(colSums(sfs), sep = "\t", file = file, append = TRUE)
  cat("\n", file = file, append = TRUE)
}

write.grapes.unfolded <- function(genes, file, nbChr = Number_of_Isolates, prefix = input_folder, title, dataset) {
  sfs <- get.sfs.grapes.unfolded(genes, nbChr, prefix)
  write.sfs.grapes.unfolded(sfs, file, nbChr, title, dataset)
}

# load list of genes in the output folder from popstat
list_of_files_POPSTAT <- list.files(path = input_folder, pattern = "*.codon.csv")
geneset_POPSTAT <- data.frame(Name = substring(list_of_files_POPSTAT, 1, nchar(list_of_files_POPSTAT) - 10), stringsAsFactors = FALSE)

# load information about effector
effector_list <- read.table(Effector_list_file, stringsAsFactors = FALSE)$V1
effector_list <- paste0("NT_",effector_list, sep="") # only if names need to be modified to match

# filter for effectors
NonEffector_POPSTAT <- subset(geneset_POPSTAT, ! Name %in% effector_list)
OnlyEffector_POPSTAT <- subset(geneset_POPSTAT, Name %in% effector_list)

# get number of genes
Total_number_genes<- as.integer(nrow(geneset_POPSTAT))
Number_NonEffectors<- as.integer(nrow(NonEffector_POPSTAT))
Number_Effectors<- as.integer(nrow(OnlyEffector_POPSTAT))

File_total_genes=paste0(Species_is,".",Total_number_genes,"_allGenes.sfs")
File_NonEffectors=paste0(Species_is,".",Number_NonEffectors,"_NonEffectors.sfs")
File_Effectors=paste0(Species_is,".",Number_Effectors,"_Effectors.sfs")

# print checks
print(paste0("The species is ",Species_is))
print(paste0("There are ",Number_of_Isolates," isolates"))
print(paste0("There are ", Total_number_genes, " genes, from which ",Number_NonEffectors, " are NonEffectors and ",Number_Effectors, " are effectors"))


#
# save sfs into a file. This files will generate the first grapes run to determine the best model per category.
write.grapes.unfolded(geneset_POPSTAT$Name,
                      file = File_total_genes,
                      title = Species_is,
                      dataset = File_total_genes)

write.grapes.unfolded(NonEffector_POPSTAT$Name,
                      file = File_NonEffectors,
                      title = Species_is,
                      dataset = File_NonEffectors)

write.grapes.unfolded(OnlyEffector_POPSTAT$Name,
                      file = File_Effectors,
                      title = Species_is,
                      dataset = File_Effectors)



##########
# part 2
##########

# We conduct bootstraps to get confidence intervals of estimates (grapes run after best model is known)
sfs.zt.all <- get.sfs.grapes.unfolded(geneset_POPSTAT$Name, nbChr = Number_of_Isolates, prefix = input_folder)
sfs.zt.NonEffector <- get.sfs.grapes.unfolded(NonEffector_POPSTAT$Name, nbChr = Number_of_Isolates, prefix = input_folder)
sfs.zt.OnlyEffector <- get.sfs.grapes.unfolded(OnlyEffector_POPSTAT$Name, nbChr = Number_of_Isolates, prefix = input_folder)


#dir.create("1_Bootstrap")

#### 
#### function for bootstrap
#### 

bootstrap <- function(sfs, nboots, prefix) {
  for (i in 1:nboots) {
    sfs.rep <- sfs[sample.int(nrow(sfs), replace = TRUE), ]
    write.sfs.grapes.unfolded(sfs.rep, nbChr = Number_of_Isolates,
                              paste("../2_bootstrap/", prefix, "_rep", i, ".sfs", sep = ""),
                              title = Species_is,
                              dataset = paste(prefix, "_rep", i, sep = ""))
  }
}

File_total_genes=paste0(Species_is,".",Total_number_genes,"_allGenes")
File_NonEffectors=paste0(Species_is,".",Number_NonEffectors,"_NonEffectors")
File_Effectors=paste0(Species_is,".",Number_Effectors,"_Effectors")

bootstrap(sfs.zt.all, 100, prefix = File_total_genes)
bootstrap(sfs.zt.NonEffector, 100, prefix = File_NonEffectors)
bootstrap(sfs.zt.OnlyEffector, 100, prefix = File_Effectors)


##########
# part 3
##########
#dir.create("2_Permutations")

geneset_POPSTAT$IsEffector <- geneset_POPSTAT$Name %in% effector_list
table(geneset_POPSTAT$IsEffector) # 117 effectors

categories <- character(length = nrow(sfs.zt.all))
#categories[geneset_POPSTAT$Name] <- "allgenes"
categories[geneset_POPSTAT$IsEffector == "TRUE"] <- "effector"
categories[geneset_POPSTAT$IsEffector == "FALSE"] <- "notEffector"
table(categories) # 117 effector and 5536 notEffector


#### 
#### function for Permutations
#### 

permute <- function(sfs, nrep, categories, prefix) {
  pb <- txtProgressBar(0, nrep, style = 3)
  for (i in 1:nrep) {
    setTxtProgressBar(pb, i)
    sfs.rep <- sfs[sample.int(nrow(sfs), replace = FALSE), ]
    for (j in unique(categories)) {
      sfs.rep.cat <- sfs.rep[which(categories == j), ]
      write.sfs.grapes.unfolded(sfs.rep.cat, nbChr = Number_of_Isolates,
                                paste("../3_permutation/", prefix, "_rep", i, "_", j, ".sfs", sep = ""),
                                title = Species_is,
                                dataset = paste(prefix, "_rep", i, "_", j, sep = ""))
    }
  }
}

# run it
Species_permutation=paste0(Species_is,"_unfolded_permutation")
permute(sfs.zt.all, 1000, categories, prefix = Species_permutation)






