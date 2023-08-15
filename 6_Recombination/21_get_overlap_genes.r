################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
######################################################################################
#   Script to take rho intervals and average in 20kb windows as in Grandaubert 2019  #
#                                                                                    #
# Files needed:                                                                      #
#             * Output files from LDhat parsed in the previous script                #
#                                                                                    #
######################################################################################
#
#
#
#####################################################################
# 


Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=105

## fixed

# create directory
working_Dir_path=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/0_scripts/15_LDhat/2_LDhat_genes")
dir.create(working_Dir_path)
setwd(working_Dir_path)

# extract information from the LDhat merged and summarized output
copy_gff=paste0("cp /home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/4_maffilter/0_files_maffilter/*.checked.gff .")# copy gff
system(copy_gff)

####################################################################################################################
#                                                   !!!IMPORTANT!!!                                                #
#                                                                                                                  #
# >> THIS PART MIGHT CHANGE ACCORDING TO THE GFF IN THE SPECIES. IT IS ALREADY DIFFERENT BTW ZYMO AND STAGO GFF << #
#                                                                                                                  #
#                                                                                                                  #
####################################################################################################################

Get_bedFile_gene=paste0("grep 'transcript' *.checked.gff | cut -d$'\t' -f1,4,5,7,8 > PerGene.txt") # get bed file for each gene in the gff
system(Get_bedFile_gene) # check PerGene.txt file !!!!!!!!!!

####################################################################################################################

# load lib
library(dplyr)

# The PerGene.txt file is generated from the specie's GFF file using bash.
# load per gene information
per_gene_df <- read.table("PerGene.txt", header = FALSE, stringsAsFactors = FALSE)

# Remove all before and up to ":".
#per_gene_df$V5 <- gsub(".*mrna.JI435_","",per_gene_df$V5)
#per_gene_df$V5 <- gsub("\\;.*","",per_gene_df$V5)
per_gene_df$V5 <- gsub("ID=","",per_gene_df$V5)

# rename columns
names(per_gene_df)[1] <- "chr"
names(per_gene_df)[2] <- "start"
names(per_gene_df)[3] <- "end"
names(per_gene_df)[4] <- "strand"
names(per_gene_df)[5] <- "geneID"

# load LDhat output
LDhat_output=paste0("/home/pereira/2020_POSTDOC_MAIN/",Species_is,"/3_analysis/13_Ldhat/7_summ/Final.",Species_is,".LDhat.20kb.txt")
chr_rho_raw <- read.table(LDhat_output, header = TRUE, stringsAsFactors = FALSE)
#chr_rho_raw$rho_mean <- na_if(chr_rho_raw$rho_mean, 0) # 0 were added and should be changed to NA since LDhat did not have a value for these windows

# create start and end columns
library(stringr)
#str_split_fixed(chr_rho_raw$window, "-", 2)[,1]
chr_rho_raw$start <- str_split_fixed(chr_rho_raw$window, "-", 2)[,1]
chr_rho_raw$end <- str_split_fixed(chr_rho_raw$window, "-", 2)[,2]

# so it is also 0-based counting as in the gff
chr_rho_raw$start <- as.integer(chr_rho_raw$start)
chr_rho_raw$end <- as.integer(chr_rho_raw$end)
chr_rho_raw$start <- chr_rho_raw$start + 1
chr_rho_raw$end <- chr_rho_raw$end + 1

# rename and order dataframes
#chr_rho_GR <- chr_rho_raw
#chr_rho_GR$window <- NULL
#col_order <- c("contig", "start", "end","rho_mean")
#chr_rho_GR <- chr_rho_GR[, col_order]
names(chr_rho_raw)[2] <- "chr"


# Compute mean recombination rates in for each gene by averaging over overlapping windows:
library(GenomicRanges)
g.rec <- makeGRangesFromDataFrame(chr_rho_raw, ignore.strand = TRUE, keep.extra.columns = TRUE)
g.gen <- makeGRangesFromDataFrame(per_gene_df, ignore.strand = TRUE)

pb <- txtProgressBar(1, nrow(per_gene_df), style = 3)
for (i in 1:nrow(per_gene_df)) {
  setTxtProgressBar(pb, i)
  r <- pintersect(g.rec, g.gen[i], drop.nohit.ranges = TRUE)
  per_gene_df[i, "rho_mean"] <- weighted.mean(r$rho_mean, r@ranges@width)
}

# export
write.table(per_gene_df,"Rho_summari_PerGene.txt", quote=F, sep="\t", row.names = F, col.names=T)


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

