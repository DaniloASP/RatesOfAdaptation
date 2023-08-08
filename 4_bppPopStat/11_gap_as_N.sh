########################################################################################
#### This script is intended to copy final set of genes, replace gapes as - to N   #####
########################################################################################

Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/1_loose_kappa

### no modification needed
#
input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/3_filtering"
#input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/Candida_auris/3_analysis/0_scripts/13_popstat/1_loose_kappa/test"
output_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/4_final_gene_set"

##############################################
### Part 1 - Copy files from MACSE NT run1 ###
##############################################

echo -e '#!/bin/bash' > copy_run1.sh
echo "input_FOLDER=$input_FOLDER" >> copy_run1.sh
echo "output_FOLDER=$output_FOLDER" >> copy_run1.sh

echo "cp $input_FOLDER/*.fasta $output_FOLDER" >> copy_run1.sh

####################################################
### Part 2 - get output sequence from each gene  ###
####################################################
# A key file is present locally, generated from the python threshold script, the file was made from run 1 aPB run 2 aPB is named the same `4_file_outgroup_KEY.txt`. Merge it aPB get unique values
#cat 4_* > 0_PN_gene_outgroup.txt
#sort -u 0_PN_gene_outgroup.txt > 0_genes_popstat_run1.txt # it has 4256, 2 more then in 7_geneset_final (has 4254), because I cancelled to jobs running too long. Check entry on # 06 JUL 2021

# do not use to one above, use this one:

echo 'for fasta_file in $output_FOLDER/*fasta' >> copy_run1.sh 
echo 'do' >> copy_run1.sh
echo 'basename $fasta_file >> 1_gene_name.txt' >> copy_run1.sh
echo 'grep '"'>g'"' $fasta_file >> 2_outgroup.txt' >> copy_run1.sh
echo 'done' >> copy_run1.sh

# prepare file
echo 'cat 1_gene_name.txt | cut -d'"'/'"' -f9 > 1_gene_nameONLY.txt' >> copy_run1.sh
echo 'sed -i '"'s/.fasta//g'"' 1_gene_nameONLY.txt' >> copy_run1.sh
echo 'sed -i '"'s/NT_//g'"' 1_gene_nameONLY.txt' >> copy_run1.sh
echo 'sed -i '"'s/>//g'"' 2_outgroup.txt' >> copy_run1.sh

echo 'paste -d "|" 1_gene_nameONLY.txt 2_outgroup.txt > popstat_gene_outgroup.txt' >> copy_run1.sh

echo "cp popstat_gene_outgroup.txt /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/" >> copy_run1.sh

##############
# Part 4     #
##############
# bpppopstat will not run with - as gaps, thus replace by N for each fasta.

# After MACSE it is necessary to replace "-" by "N" using `sed -i '' 's/-/N/g' *.fasta`. `grep -l` prints one hit per file.

echo 'for file in $output_FOLDER/*fasta' >> copy_run1.sh
echo 'do' >> copy_run1.sh
echo 'sed -i '"'s/-/N/g'"' $file' >> copy_run1.sh
echo 'done' >> copy_run1.sh


# wallace
sbatch --job-name=VD_sct11 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=NONE --mail-user=pereira@evolbio.mpg.de --partition=fast < copy_run1.sh

