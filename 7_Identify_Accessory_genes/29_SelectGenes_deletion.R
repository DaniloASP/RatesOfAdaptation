################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
######################################################################################
#   This script will overlap gene coordinates with deletion coordinates              #
#                                                                                    #
# Files needed:                                                                      #
#             * Output from CNVnator                                                 #
#                                                                                    #
######################################################################################
#
#
#
#####################################################################
# Analyse output from CNVnator from various species

Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=105


# set wdir
setwd("~/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp")

### summarize CNVnator output
# set path to CNVnator output and grep all files
grep_raw <- paste0("grep '' /home/pereira/2020_POSTDOC_MAIN/",{Species_is},"/3_analysis/15_CNVnator/2_output/*.CNV.output > ",{Abb_slurm},"_CNV.raw")
system(grep_raw, intern = TRUE)
###

# lib load
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(ape)
library(GenomicFeatures)
#BiocManager::install("Repitools")
library(Repitools)
library(BiocGenerics)
library(BiocManager)
library(parallel)
library(data.table)


######## load GFF
# extract information from the LDhat merged and summarized output
find_gff=paste0("find /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/4_maffilter/0_files_maffilter/*.checked.gff")
GFF_species <- system(find_gff, intern = TRUE)

# load GFF file or bed file with chr / start / end / gene_name
reference_GFF <- read.gff(GFF_species, na.strings = c(".", "?"), GFF3 = TRUE)

# gene names have to match those in
Popstat_Output=paste0("ls /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/15_popstat/3_fixed_kappa/*log")
system(Popstat_Output)

# work on gff so gene name match popstat output
table(reference_GFF$type)

gene_name_only <- reference_GFF %>% filter(type == "mRNA") # OBS: this will change per species
table(reference_GFF$attributes)[1]
gene_name_only$attributes <- gsub(";.*","",gene_name_only$attributes)
gene_name_only$attributes <- gsub("ID=rna-","",gene_name_only$attributes)

#gene_name_only$attributes <- gsub("ID=rna-","",gene_name_only$attributes)
#gene_name_only$attributes <- gsub("RCO","RC0",gene_name_only$attributes)

#gene_name_only$attributes <- paste(gene_name_only$attributes,"A", sep = "")
#gene_name_only$seqid <- reorder(gene_name_only$seqid)
table(gene_name_only$seqid)
gene_name_only$attributes

# if names in the popstat output folder and GFF match, create a genomicrange GFF object
gr_gff <- GRanges(
  seqnames = Rle(gene_name_only$seqid), #chromosome data here
  ranges = IRanges(gene_name_only$start, end = gene_name_only$end, names = (gene_name_only$attributes)), #size range and id
  strand = Rle((gene_name_only$strand)), # strand orientation + or -
  gene = gene_name_only$attributes) # isolate id


######## load raw CNVnator output
# load raw greped file and simplify first collumn
print(paste0({Abb_slurm},"_CNV.raw"))
CNVnator_DF <- read.table(paste0({Abb_slurm},"_CNV.raw"))
CNVnator_DF$V1 <- gsub(".*2_output/","",CNVnator_DF$V1) # remove 2_output/ plus everything before it

# split first column for isolates
CNVnator_DF <- CNVnator_DF %>%
  separate(V1, c("isolate", "CNV_type"), ".CNV.output:")

# CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0
colnames(CNVnator_DF) <- c("isolate", "CNV_type", "coordinates","CNV_size", "normalized_RD","e_val1", "e_val2","e_val3", "e_val4","q0")

# split first column to make a new column with coordinates
CNVnator_DF$coordinates_raw <- CNVnator_DF$coordinates

CNVnator_DF <- CNVnator_DF %>%
  separate(coordinates_raw, c("chromosome", "range"), ":")

CNVnator_DF <- CNVnator_DF %>%
  separate(range, c("start", "end"), "-")

# Filters used as in Hartmann 2017 and Stalder 2021. OBS: Hartmann 2017 validated gene deletions with PCR
PASS_CNV_DF <- CNVnator_DF %>% filter(e_val1 < 0.05 & q0 < 0.5) %>% filter(CNV_type == "deletion" & normalized_RD < 0.4 | CNV_type == "duplication" & normalized_RD > 1.6)

table(PASS_CNV_DF$CNV_type)

# save PAV numbers
write.table(table(PASS_CNV_DF$CNV_type), paste0({Abb_slurm},".numberPAV.txt"), quote = FALSE, row.names = FALSE)

# select only deletions
PASS_deletions_DF <- PASS_CNV_DF %>% filter(CNV_type == "deletion")
PASS_deletions_DF$start = as.numeric(PASS_deletions_DF$start)
PASS_deletions_DF$end = as.numeric(PASS_deletions_DF$end)

# Create a genome range object to use a intersection function later. Note that there should be no overlap in deletions for a same isolate
gr_deletion <- GRanges(
  seqnames = Rle(PASS_deletions_DF$chromosome), #chromosome data here
  ranges = IRanges(PASS_deletions_DF$start, end = PASS_deletions_DF$end, names = (PASS_deletions_DF$coordinates)), #size range and deletion id
  isolate = PASS_deletions_DF$isolat) # isolate id

# Select genes completly inside deleted regions
# intersect
#findOverlaps(gr_gff,gr_deletion, type = "within")
#intersect(gr_gff,gr_deletion)
overlap_gr <- mergeByOverlaps(gr_gff, gr_deletion, type = "within") # https://www.biostars.org/p/9474374/

# slice to export
export_df <- overlap_gr[,2]
export_df <- as.data.frame(export_df)
colnames(export_df)[1] <- "gene"
export_df$isolate <- overlap_gr[,4]
#write.table(export_df, "Gene_deletion_ZT_typeWithin.txt", quote = FALSE, row.names = FALSE)
overlap_gr_df<-data.frame(lapply(overlap_gr, as.character), stringsAsFactors=FALSE) # https://stackoverflow.com/questions/48709160/why-do-i-keep-getting-an-error-when-creating-a-csv-file-in-r
#write.table(data1, "Full_Gene_deletion_ZT_typeWithin.txt", quote = FALSE, row.names = FALSE)

# check distribution of deletions.
gene_deletion_count <- as.data.frame(table(export_df$gene))
gene_deletion_count$del_perc <- (gene_deletion_count$Freq / Num_chr) * 100
hist(gene_deletion_count$del_perc)

# from GFF annotation, remove genes showing impartial deletion in the dataset
CORE_gene_set_DF <- gene_name_only[!(gene_name_only$attributes %in% gene_deletion_count$Var1),]
table(CORE_gene_set_DF$seqid)

genes_in_grapes_run <- system(Popstat_Output, intern = TRUE)
genes_in_grapes_run
CORE_gene_set_DF$attributes
genes_in_grapes_run <- gsub(".*NT_", "",genes_in_grapes_run)
genes_in_grapes_run <- gsub(".log", "",genes_in_grapes_run)
#genes_in_grapes_run <- gsub("RCO7", "RC07_",genes_in_grapes_run)
genes_in_grapes_run

# From list of core genes, select only the core used in the grapes run
CORE_grapes_set_DF <- CORE_gene_set_DF[(CORE_gene_set_DF$attributes %in% genes_in_grapes_run),]
table(CORE_grapes_set_DF$seqid)


# Export file with genes in deletion region
write.table(gene_deletion_count, paste0({Abb_slurm},".GenePAV.txt"), quote = FALSE, row.names = FALSE)
write.table(export_df, paste0({Abb_slurm},".GeneAndIsolate_PAV.txt"), quote = FALSE, row.names = FALSE)


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

