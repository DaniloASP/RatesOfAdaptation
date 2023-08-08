# Functional analysis of proteins 

################
# interproscan #
################
Species_is="Verticillium"
Abb_is="Vdahliae"


# copy data
ls /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan/
cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan

cp /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/8_select_orthologs/*_protein.fasta /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan

# first remove asterix from the protein file or it will not run (was necessary only for ND)
#sed -i 's/*//g' /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan/Asp_NRRL3357_protein.simpleHeader.fasta

# run.
module load java/x64/11u1
module load perl/5.30.1
module load python/3.7.1
echo -e '#!/bin/bash' > interpro.run.sh
echo "/data/biosoftware/InterProScan/interproscan-5.48-83.0/interproscan.sh -appl SignalP_EUK-4.1 -appl TMHMM-2.0c -dp --goterms --pathways -iprlookup -f TSV, GFF3 -i Vdahliae_protein.fasta" >> interpro.run.sh

# wallace
sbatch --job-name=VD_sc8 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=64G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=standard < interpro.run.sh

############################################# continue here

Species_is="Verticillium"
Abb_is="Vdahliae"

output_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/3_fixed_kappa"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan

#########
module load python/python-anaconda-5.0.1 
source activate Danilo_postdoc

module load R/4.1.2

nano 1_functional_ann.r # execute Rscript 1_functional_ann.r

#################################################################### copy and paste in 1_functional_ann.r
#####Functional annotation########

# PART 1
# THIS FIRST PART IS TO GET THE OUTPUT FROM INTERPROSCAN AND RESHAPE IT INTO INPUT FOR GO ENRICHMENT TEST

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

#################################################################### paste end
# execute Rscript 1_functional_ann.r

# select proteins found in secretome
cut -d',' -f1 functional_anno_secretome.csv > functional_anno_secretome.names.csv # open and remove "

mv functional_anno_secretome.names.csv secretome.names.tmp.csv && tail -n +2 secretome.names.tmp.csv > functional_anno_secretome.names.csv

sed -i'' 's/"//' functional_anno_secretome.names.csv && sed -i'' 's/"//' functional_anno_secretome.names.csv

seqtk subseq *_protein.fasta functional_anno_secretome.names.csv > name.protein.secretome.fa

# remove whitespace from fasta header and create a local copy
sed 's/ .*//' name.protein.secretome.fa > name.protein.secretome.fixed.fa

python /home/pereira/software/EffectorP-3.0/EffectorP.py -i name.protein.secretome.fixed.fa > name.protein.secretome.effectorP_OUT.txt


# extract effector sequences
grep 'XP_' name.protein.secretome.effectorP_OUT.txt > name.protein.secretome.effectorP_OUT.forR.txt

nano 2_effectorP.r # execute Rscript 2_effectorP.r

#################################################################### copy and paste in 2_effectorP.r
#####Functional annotation########

# PART 1
# THIS FIRST PART IS TO GET THE OUTPUT FROM EFFECTORP AND RESHAPE INTO A DATA FRAME

library(reshape)
library(dplyr)

# read in InterProScan TSV file
file <- "name.protein.secretome.effectorP_OUT.forR.txt"

effectorP_OUT <- read.csv(file, sep="\t", fill=TRUE, header=FALSE)
effectorP_OUT$V1 <- sub(" .*", "", effectorP_OUT$V1) # remove strings after space

# subset for genes with signal pep and no transmembrane helices
effectors_only_DF <- subset(effectorP_OUT, V5 != "Non-effector")

# write output
write.table(effectors_only_DF, ("name.protein.secretome.effectorP_OUT.forR.effectorONLY.txt"), sep=",",  col.names=FALSE,row.names=FALSE, na="", quote = FALSE)

#################################################################### paste end
# execute Rscript 2_effectorP.r

# get effectors name 
cut -d',' -f1 name.protein.secretome.effectorP_OUT.forR.effectorONLY.txt > effectors.namesONLY.txt

# Based on the final gene dataset, now find the effectors that survided. Easier in R.
ls $output_FOLDER/*log > 0_list_final_geneset.txt

mv 0_list_final_geneset.txt 0_list_final_geneset.temp.txt && cut -d'/' -f9 0_list_final_geneset.temp.txt > 0_list_final_geneset.txt && rm 0_list_final_geneset.temp.txt

nano 3_findEffectors.r # execute Rscript 3_findEffectors.r

#################################################################### copy and paste in 3_findEffectors.r
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
#################################################################### paste ends
# execute Rscript 3_findEffectors.r



##### If XM not matching XP for protein nomenclature, this will be necessary:
## FOR SM ONLY, because genes are coded as XM_ instead of XP_
cp /home/pereira/2020_POSTDOC_MAIN/Sphaerulina/3_analysis/2_reference_files/gff/Script10_fix.txt .

IFS=$'|'
while read line
do
linearray=( $line )
new_name=${linearray[0]}
old_name=${linearray[1]}
echo $old_name
#echo ">Ztritici:$CHR_id:1:+:$length_id"
sed -i "s/$old_name/$new_name/g" 0_list_final_geneset.txt
done < Script10_fix.txt
####

###
cp effectors.namesONLY.txt effectors.names2ONLY.txt

IFS=$'|'
while read line
do
linearray=( $line )
new_name=${linearray[0]}
old_name=${linearray[1]}
echo $old_name
#echo ">Ztritici:$CHR_id:1:+:$length_id"
sed -i "s/$old_name/$new_name/g" effectors.names2ONLY.txt
done < Script10_fix.txt

## mac
IFS=$'|'
while read line
do
linearray=( $line )
new_name=${linearray[0]}
old_name=${linearray[1]}
echo $old_name
#echo ">Ztritici:$CHR_id:1:+:$length_id"
sed -i '' -e "s/$new_name/$old_name/g" effectors.names2ONLY.txt
done < Script10_fix.txt


nano 3_findEffectors.r # execute Rscript 3_findEffectors.r
