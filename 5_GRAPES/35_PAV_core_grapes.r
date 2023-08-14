# # create categories for GRAPES run




Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=60 # number of isolates in the downsampled dataset to compare PAV and core genes
Best_grapes_model_allgenes="GammaZero"



# fixed

# load lib
library(dplyr)
library(Hmisc)
library(tidyverse)

setwd(paste0("/Users/danilo/Dropbox/Backup/MacBookPro/PostDoc/Stukenbrock/Project/2020_postdoc_main/03_results/Common_analisys/20220726_PAV_input_grapes/",Species_is,"/3_grapes"))

# create dir for output
Popstat_Output=paste0("../2_popstat/2_merged/")


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


###################################################################################
# part 2 - load data about secreted proteins and effectors and create 3 categories
###################################################################################

# load data
genesCORE <- read.table("popstat_CORE_outgroup.txt", sep="|", stringsAsFactors = FALSE, header = FALSE)
genesPAV <- read.table("popstat_PAV_outgroup.txt", sep="|", stringsAsFactors = FALSE, header = FALSE)


### thin list to the gene dataset in the species
# load list of genes in the output folder from popstat
list_of_files_POPSTAT <- list.files(path = Popstat_Output, pattern = "*.codon.csv")

# convert to dataframe
geneset_POPSTAT <- data.frame(Name = substring(list_of_files_POPSTAT, 1, nchar(list_of_files_POPSTAT) - 10), stringsAsFactors = FALSE)

# match gene names
genesCORE$V1 <- paste0("NT_", genesCORE$V1)
genesPAV$V1 <- paste0("NT_", genesPAV$V1)

# I need to remove some genes in the PAV dataset to make the number of genes equal between datasets.
remove_genesPAV <- head(genesPAV, - 3102)
geneset_POPSTAT <- geneset_POPSTAT %>% filter(!Name %in% remove_genesPAV$V1)

# set categories
geneset_POPSTAT$Category <- "CORE"

geneset_POPSTAT$Category[geneset_POPSTAT$Name %in% genesPAV$V1] <- "ACC"
table(geneset_POPSTAT$Category) 

All_categories <- unique(geneset_POPSTAT$Category)

geneset_POPSTAT$Species <- Species_is

# export
write.table(geneset_POPSTAT,paste0(Abb_slurm,"_Gene_category_list.txt"), quote=F, sep="\t", row.names = F, col.names=F)


############################################################################
# part 3 - CREATE OUTPUT FILES FOR BOOTSTRAP
############################################################################

# create a loop to extract the gene sample from each recombination category
# We conduct bootstraps to get confidence intervals of estimates (grapes run after best model is known)

# this loop will extract categories to individual dataframes, and will automatically name the dataframes according to the category
vector_of_GeneCategory <- vector()
for (GeneCategory in All_categories) {
  #print(vector_of_GeneCategory)
  name_DF <- paste("DF", GeneCategory, sep="_") # this will construct a name from within the loop
  vector_of_GeneCategory <- c(vector_of_GeneCategory, name_DF)
  #assign(name_DF, loop_DF) # this will assign to that name a dataframe
  assign(name_DF, (geneset_POPSTAT[with(geneset_POPSTAT, Category %in% GeneCategory),]))
}


# this loop will create the SFS for each gene category
for (name_DF in vector_of_GeneCategory){
  cur.df <- get(name_DF)
  assign(name_DF, get.sfs.grapes.unfolded(cur.df$Name, nbChr = Num_chr, prefix = Popstat_Output))
}

#### 
#### function for bootstrap
#### 

dir.create("0_FullDataset_Bootstrap")
bootstrap <- function(sfs, nboots, prefix) {
  pb <- txtProgressBar(0, nboots, style = 3)
  for (i in 1:nboots) {
    setTxtProgressBar(pb, i)
    sfs.rep <- sfs[sample.int(nrow(sfs), replace = TRUE), ]
    write.sfs.grapes.unfolded(sfs.rep, nbChr = Num_chr,
                              paste("./0_FullDataset_Bootstrap/", prefix, "_rep", i, ".sfs", sep = ""),
                              title = Species_is,
                              dataset = paste(prefix, "_rep", i, sep = ""))
  }
}


# this loop will create the bootstrap replicates for SFS of each gene category
n = length(vector_of_GeneCategory)
count_bar = 1
pb <- txtProgressBar(1, n, style = 3)
for (name_DF in vector_of_GeneCategory){
  cur.df <- get(name_DF)
  bootstrap(cur.df, 100, prefix = name_DF)
}

Tar_gz_file=paste0("tar -zcvf 0_FullDataset_Bootstrap.tar.gz 0_FullDataset_Bootstrap")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!


############################################################################
# part 4 - CREATE OUTPUT FILES FOR PERMUTATION
############################################################################
dir.create("1_FullDataset_Permutation")

sfs.FullDataset <- get.sfs.grapes.unfolded(geneset_POPSTAT$Name, nbChr = Num_chr, prefix = Popstat_Output)
Gene_category_perm <- geneset_POPSTAT$Category
table(Gene_category_perm)


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
permute(sfs.FullDataset, 1000, Gene_category_perm, prefix = "FullData")

Tar_gz_file=paste0("tar -zcvf 1_FullDataset_Permutation.tar.gz 1_FullDataset_Permutation")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!

############################################################################
# part 4 - CREATE OUTPUT FILES FOR SINGLE RUN OF EACH CATEGORY 
############################################################################

## NonSecreted
CORE_POPSTAT <- geneset_POPSTAT %>%
  filter(str_detect(Category, "CORE"))

## secreted
ACC_POPSTAT <- geneset_POPSTAT %>%
  filter(str_detect(Category, "ACC"))

write.grapes.unfolded(CORE_POPSTAT$Name,
                      file = "CORE.genes.sfs",
                      title = Species_is,
                      dataset = "CORE.genes")

write.grapes.unfolded(ACC_POPSTAT$Name,
                      file = "ACC.genes.sfs",
                      title = Species_is,
                      dataset = "ACC.genes")

Tar_gz_file=paste0("tar -zcvf 2_FullDataset_indv_cat.tar.gz *genes.sfs")
system(Tar_gz_file) # check PerGene.txt file !!!!!!!!!!

### end


