################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
###############################################################################
#   This script will prepare files and sequences for input in bppPopStat      #
#                                                                             #
# Files needed:                                                               #
#             * Output from script 10_Gap_threshold.v2.sh                     #
#                                                                             #
#                                                                             #
###############################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################
Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/1_loose_kappa

# Path to variables
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

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
