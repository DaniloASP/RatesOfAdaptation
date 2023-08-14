# Script to take rho intervals and average in 20kb windows as in Grandaubert 2019
# # create categories for GRAPES run

Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=105



# fixed

# load lib
library(dplyr)
library(Hmisc)

mk_nw_dir=paste0("mkdir /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/20_grapes_recombination")
system(mk_nw_dir)

working_Dir_path=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/20_grapes_recombination")
setwd(working_Dir_path)

copy_files_to_wd=paste0("cp /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/0_scripts/15_LDhat/2_LDhat_genes/Rho_summari_PerGene.txt .")
system(copy_files_to_wd)


copy_files_to_wd=paste0("cp /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/18_grapes_3genecat/*anno_secretome.names.csv .")
system(copy_files_to_wd)

copy_files_to_wd=paste0("cp /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/18_grapes_3genecat/*namesONLY.txt .")
system(copy_files_to_wd)

Popstat_Output=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/15_popstat/3_fixed_kappa/")

get_PAV_genes=paste0("cp /home/pereira/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp/",{Abb_slurm},".GenePAV.txt .")
system(get_PAV_genes)

############################################################
# part 1 - load functions
############################################################

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

write.grapes.unfolded <- function(genes, file, nbChr = Num_chr, prefix = Popstat_Output, title, dataset) {
  sfs <- get.sfs.grapes.unfolded(genes, nbChr, prefix)
  write.sfs.grapes.unfolded(sfs, file, nbChr, title, dataset)
}


###################################################################
# part 2 - load data - FULL GENE SET (NOT REMOVING EFFECTORS HERE)
###################################################################

# load data
dat <- read.table("Rho_summari_PerGene.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
Effector_file <- system("ls *namesONLY.txt", intern=T) # The effector file needs to be transfered to wallace manually
Effector_list <- read.table(Effector_file, sep="\t", stringsAsFactors = FALSE, header = FALSE)
secretome_file <- system("ls *anno_secretome.names.csv", intern = TRUE)
secretome_df <- read.table(secretome_file, sep=",", stringsAsFactors = FALSE, header = FALSE)
genesPAV_file <- paste0({Abb_slurm},".GenePAV.txt")
genesPAV_df <- read.table(genesPAV_file, stringsAsFactors = FALSE, header = TRUE)


# select only core chr
#dat_core <- dat %>% filter(chr < 14)
dat_core <- dat


#chr_rho_raw$rho_mean <- na_if(chr_rho_raw$rho_mean, 0) # 0 were added and should be changed to NA since LDhat did not have a value for these windows
dat_core$rho_mean[is.na(dat_core$rho_mean)] <- 0

# match gene names
dat_core$geneID <- paste0("NT_", dat_core$geneID) # add NT_ to begin of gene name
Effector_list$V1 <- paste0("NT_", Effector_list$V1) # add NT_ to begin of gene name
secretome_df$V1 <- paste0("NT_", secretome_df$V1)
genesPAV_df$Var1 <- paste0("NT_", genesPAV_df$Var1)

### thin list to the gene dataset in the species
# load list of genes in the output folder from popstat
list_of_files_POPSTAT <- list.files(path = Popstat_Output, pattern = "*.codon.csv")

# convert to dataframe
geneset_POPSTAT <- data.frame(Name = substring(list_of_files_POPSTAT, 1, nchar(list_of_files_POPSTAT) - 10), stringsAsFactors = FALSE)
#dat_core$geneID <- gsub("NT_RCO7_", "NT_RCO7", dat_core$geneID) # the name should match geneset_POPSTAT$Name

# thin dataset
# if gene is in the list of PAV, it will be removed
geneset_POPSTAT <- geneset_POPSTAT %>% filter(!Name %in% genesPAV_df$Var1) # will remove genes showing PAV, thus keeping only core genes
geneset_POPSTAT <- geneset_POPSTAT %>% filter(!Name %in% secretome_df$V1) # will remove genes classified as secreted
geneset_POPSTAT <- geneset_POPSTAT %>% filter(!Name %in% Effector_list$V1) # will remove genes classified as effectors

dat_core_final <- subset(dat_core, geneID %in% geneset_POPSTAT$Name)

# categorize

dat_core_final$RHO_CLASS <- cut2(dat_core_final$rho_mean, g = 2, levels.mean = TRUE) # note g=20 here, but there were created 15 categories
#rec <- as.numeric(levels(dat_core$RHO_CLASS))
rec <- unique(dat_core_final$RHO_CLASS)
genes_per_rho <- as.data.frame(table(dat_core_final$RHO_CLASS))
names(genes_per_rho)[1] <- "Rho"
names(genes_per_rho)[2] <- "N_of_genes_total"

table(dat_core_final$RHO_CLASS)


write.table(dat_core_final, paste0(Abb_slurm,"-Core_noSecr_noEff_genePerGroup_rho15.txt"), quote = FALSE, row.names = FALSE)   # export
write.table(genes_per_rho, paste0(Abb_slurm,"-Core_noSecr_noEff_Freq_rho15.txt"), quote = FALSE, row.names = FALSE)   # export

#####################################
# part 3 - CREATE SFS FOR BOOTSTRAP #
#####################################

# create a loop to extract the gene sample from each recombination category
# We conduct bootstraps to get confidence intervals of estimates (grapes run after best model is known)

# this loop will extract categories to individual dataframes, and will automatically name the dataframes according to the category
vector_of_DF <- vector()
for (rec_category in rec) {
  #print(rec_category)
  name_DF <- paste("DF", rec_category, sep="_") # this will construct a name from within the loop
  vector_of_DF <- c(vector_of_DF, name_DF)
  #assign(name_DF, loop_DF) # this will assign to that name a dataframe
  assign(name_DF, (dat_core_final[with(dat_core_final, RHO_CLASS %in% rec_category),]))
}


# this loop will create the SFS for each gene category
for (name_DF in vector_of_DF){
  cur.df <- get(name_DF)
  
  assign(name_DF, get.sfs.grapes.unfolded(cur.df$geneID, nbChr = Num_chr, prefix = Popstat_Output))
}

#### 
#### function for bootstrap
#### 

dir.create("0_FullDataset_Bootstrap")
bootstrap <- function(sfs, nboots, prefix) {
  for (i in 1:nboots) {
    sfs.rep <- sfs[sample.int(nrow(sfs), replace = TRUE), ]
    write.sfs.grapes.unfolded(sfs.rep, nbChr = Num_chr,
                              paste("./0_FullDataset_Bootstrap/", prefix, "_rep", i, ".sfs", sep = ""),
                              title = Species_is,
                              dataset = paste(prefix, "_rep", i, sep = ""))
  }
}


# this loop will create the bootstrap replicates for SFS of each gene category
n = length(vector_of_DF)
count_bar = 1
pb <- txtProgressBar(1, n, style = 3)
for (name_DF in vector_of_DF){
  count_bar = count_bar + 1
  setTxtProgressBar(pb, count_bar)
  cur.df <- get(name_DF)
  bootstrap(cur.df, 100, prefix = name_DF)
}

Tar_gz_file=paste0("tar -zcvf 0_FullDataset_Bootstrap.tar.gz 0_FullDataset_Bootstrap")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!


############################################################################
# part 4 - CREATE OUTPUT FILES - GENE SET WITH EFFECTORS
############################################################################
dir.create("1_FullDataset_Permutation")

sfs.FullDataset <- get.sfs.grapes.unfolded(dat_core_final$geneID, nbChr = Num_chr, prefix = Popstat_Output)
categories <- character(length = nrow(sfs.FullDataset))
categories <- as.vector(dat_core_final$RHO_CLASS)
table(categories)


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
      write.sfs.grapes.unfolded(sfs.rep.cat, nbChr = Num_chr,
                                paste("./1_FullDataset_Permutation/", prefix, "_rep", i, "_", j, ".sfs", sep = ""),
                                title = Species_is,
                                dataset = paste(prefix, "_rep", i, "_", j, sep = ""))
    }
  }
}


# run it
# this loop will create the bootstrap replicates for SFS of each gene category
permute(sfs.FullDataset, 1000, categories, prefix = "FullData")

Tar_gz_file=paste0("tar -zcvf 1_FullDataset_Permutation.tar.gz 1_FullDataset_Permutation")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!

############################################################################
# part 5 - CREATE OUTPUT FILES - FULL GENE SET (NOT REMOVING EFFECTORS HERE)
############################################################################
dir.create("2_RealData_mean")

# create a loop to extract the gene sample from each recombination category
# We conduct bootstraps to get confidence intervals of estimates (grapes run after best model is known)

# this loop will extract categories to individual dataframes, and will automatically name the dataframes according to the category
vector_of_DF <- vector()
for (rec_category in rec) {
  #print(rec_category)
  name_DF <- paste("DF", rec_category, sep="_") # this will construct a name from within the loop
  vector_of_DF <- c(vector_of_DF, name_DF)
  #assign(name_DF, loop_DF) # this will assign to that name a dataframe
  assign(name_DF, (dat_core_final[with(dat_core_final, RHO_CLASS %in% rec_category),]))
}


# this loop will create the SFS for each gene category. This is the first run necessary to find the best model per category
for (name_DF in vector_of_DF){
  cur.df <- get(name_DF)
  print(name_DF)
  file_name <- paste("./2_RealData_mean/", Abb_is, "_", name_DF, ".sfs", sep="")
  dataset_name <- paste(Abb_is, "_", name_DF, sep="")
  write.grapes.unfolded(cur.df$geneID, file = file_name, title = Abb_is, dataset = dataset_name)
}

Tar_gz_file=paste0("tar -zcvf 2_RealData_mean.tar.gz 2_RealData_mean")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!
