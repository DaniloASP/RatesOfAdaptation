################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################################
# This script is intended to filter alignments after MACSE                         #
#                                                                                  #
# Files needed:                                                                    #
#             * Folder with MACSE output alignment files                           #
#                                                                                  #
#                                                                                  #
####################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################


# 1 - Copy output from MACSE NT folder
# 2 - Get genes with frameshift (coded as !)
# 3 - Get genes with any gap at the beginning of alignment. (to be remove later) OR run macse penalizing gaps at terminal parts (favor gaps inside alignment)
# 4 - Get genes with any gap at end of alignment. (to be remove later) OR run macse penalizing gaps at terminal parts (favor gaps inside alignment)
# 5 - Get files to perform quantification of within sequence gaps.

Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/12_filtering

#########################################
### Part 1 - Copy files from MACSE NT ###
#########################################
echo -e '#!/bin/bash' > scp9_filter.sh
echo "cp /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/2_macse_output/1_NT/*.fasta /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/3_filtering" >> scp9_filter.sh

####################################
# Part 2 - Frameshift filtering    #
####################################
input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/3_filtering"

# Get genes with frame shift
echo "grep -l '"'!'"' $input_FOLDER/*fasta >> remove_genes_list.txt" >> scp9_filter.sh

######################################
### Part 3 - Terminal GAPS (start) ###
######################################
echo "grep -l '"'^-'"' $input_FOLDER/*fasta >> remove_genes_list.txt" >> scp9_filter.sh

######################################
### Part 4 - Terminal GAPS (end)   ###
######################################
echo "grep -l '"'\-$'"' $input_FOLDER/*fasta >> remove_genes_list.txt" >> scp9_filter.sh

######################################
### Part 5 - START CODON           ###
######################################
echo "grep -L '"'^atg'"' $input_FOLDER/*fasta >> remove_genes_list.txt" >> scp9_filter.sh

echo 'while read line; do rm $line; done < remove_genes_list.txt' >> scp9_filter.sh



sbatch --job-name=VD_sct9 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=NONE --mail-user=pereira@evolbio.mpg.de --partition=global scp9_filter.sh


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
