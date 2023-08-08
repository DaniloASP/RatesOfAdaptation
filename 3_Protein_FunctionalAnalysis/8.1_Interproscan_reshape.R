################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################################
# This script will annotate functions to the various protein sequences             #
#                                                                                  #
# Files needed:                                                                    #
#             * Target species protein fasta                                       #
#                                                                                  #
# Tools manuals: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html  #
#                                                                                  #
####################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

library(reshape)
library(dplyr)

# read in InterProScan TSV file
file <- system("ls *.fasta.tsv", intern=T)

interpro.raw <- read.csv(file, sep="\t", fill=TRUE, header=F, col.names=1:15)

names(interpro.raw) <- c("gene_id", "MD5", "length", "analysis", "signature_acc", "signature_desc", "signature_start", "signature_stop", "signature_score", "signature_status", "IPR_date", "IPR_accession", "IPR_description", "GO_terms", "Pathway")

### Summarize SignapP predictions
SignalPcounts.df <- as.data.frame(interpro.raw %>% group_by(gene_id) %>% filter(analysis == "SignalP_EUK") %>% summarize(SignalP = n()))

### Summarize TMHMM
TMHMMcounts.df <- as.data.frame(interpro.raw %>% group_by(gene_id) %>% filter(analysis == "TMHMM", signature_acc == "TMhelix") %>% summarize(TMHMM = n()))

### Merge all data
df.list <- list(SignalPcounts.df, TMHMMcounts.df)

match.by <- "gene_id"
annot.df = Reduce(function(...) merge(..., by=match.by, all=T), df.list)

# subset for genes with signal pep and no transmembrane helices (= NA in TMHMM and 1 in SignalP)
annot.secretome.df <- subset(annot.df, is.na(TMHMM) & SignalP == 1)

# write output
write.csv(annot.df, ("functional_anno_full.csv"), row.names=F, na="")
write.csv(annot.secretome.df, ("functional_anno_secretome.csv"), row.names=F, na="")

# select proteins found in secretome
system("cut -d',' -f1 functional_anno_secretome.csv > functional_anno_secretome.names.csv")

system("mv functional_anno_secretome.names.csv secretome.names.tmp.csv && tail -n +2 secretome.names.tmp.csv > functional_anno_secretome.names.csv")

system("sed -i'' 's/"//' functional_anno_secretome.names.csv && sed -i'' 's/"//' functional_anno_secretome.names.csv")

system("seqtk subseq *_protein.fasta functional_anno_secretome.names.csv > name.protein.secretome.fa")

# remove whitespace from fasta header and create a local copy
system("sed 's/ .*//' name.protein.secretome.fa > name.protein.secretome.fixed.fa")

system("python /home/pereira/software/EffectorP-3.0/EffectorP.py -i name.protein.secretome.fixed.fa > name.protein.secretome.effectorP_OUT.txt")

# extract effector sequences
system("grep 'XP_' name.protein.secretome.effectorP_OUT.txt > name.protein.secretome.effectorP_OUT.forR.txt")

# read in InterProScan TSV file
file <- "name.protein.secretome.effectorP_OUT.forR.txt"

effectorP_OUT <- read.csv(file, sep="\t", fill=TRUE, header=FALSE)
effectorP_OUT$V1 <- sub(" .*", "", effectorP_OUT$V1) # remove strings after space

# subset for genes with signal pep and no transmembrane helices
effectors_only_DF <- subset(effectorP_OUT, V5 != "Non-effector")

# write output
write.table(effectors_only_DF, ("name.protein.secretome.effectorP_OUT.forR.effectorONLY.txt"), sep=",",  col.names=FALSE,row.names=FALSE, na="", quote = FALSE)

# get effectors name 
system("cut -d',' -f1 name.protein.secretome.effectorP_OUT.forR.effectorONLY.txt > effectors.namesONLY.txt")

# Based on the final gene dataset, now find the effectors that survided. Easier in R.
system("ls $output_FOLDER/*log > 0_list_final_geneset.txt")

system("mv 0_list_final_geneset.txt 0_list_final_geneset.temp.txt && cut -d'/' -f9 0_list_final_geneset.temp.txt > 0_list_final_geneset.txt && rm 0_list_final_geneset.temp.txt")


# count number of effectors in the final gene dataset after popstat fixed kappa run

library(reshape)
library(dplyr)

# read in InterProScan TSV file
Genedataset <- "0_list_final_geneset.txt"
Effectorsdata<- "effectors.namesONLY.txt"

Genedataset_df <- read.csv(Genedataset, sep="\t", fill=TRUE, header=F)
Effectordataset_df <- read.csv(Effectorsdata, sep="\t", fill=TRUE, header=F)

# clean names
Genedataset_df$V1 <- gsub("NT_","", Genedataset_df$V1)
Genedataset_df$V1 <- gsub(".log","", Genedataset_df$V1)

# mark effectors
Genedataset_df$V2 <- (as.vector(Genedataset_df$V1) %in% as.vector(Effectordataset_df$V1))
as.data.frame(table(Genedataset_df$V2)) # true means match with effector list, thus TRUE = EFFECTOR. in CP 3078 FALSE/ 85 TRUE

write.table(as.data.frame(table(Genedataset_df$V2)), ("NonEffectors_effectors_number.txt"), sep=",",  col.names=FALSE,row.names=FALSE, na="", quote = FALSE)


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
