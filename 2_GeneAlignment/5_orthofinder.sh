################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
###########################################################################
# This script will prepare all necessary files for pseudogenome assembly  #
#                                                                         #
# Files needed:                                                           #
#               * protein fasta file from target and outgroup species     #
#                                                                         #
###########################################################################
#
#
#
########################################### Script will do ########################################### 
# File input needed for pseudogenome preparation                                                     #
# * Run orthofinder to identify unique ortologue gene copies between target and outgroup species     #
#                                                                                                    #
######################################################################################################              

######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

## change needed per species
Species_is="Verticillium"
Abb_is="Vdahliae"

## copy to orthofinder folder
BASE="/home/pereira/2020_POSTDOC_MAIN/$Species_is"
REF="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/0_references/0_genome_DNA/$Abb_is"

Output_folder="$BASE/3_analysis/8_select_orthologs"

cd $Output_folder

# target species
cp "$REF"_protein.fasta $Output_folder

# outgroup species
cp /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/4_BRAKER2/braker/augustus.hints.aa $Output_folder/augustus.hints.fasta

# simplify the header in the protein fasta file
sed -i'' 's/ uncharacterized.*//' "$Abb_is"_protein.fasta
sed -i'' 's/ LOW.*//' "$Abb_is"_protein.fasta

echo $(pwd)

# activate conda environment
conda activate env_orthofinder

# run it. after -f it should be a folder containing the protein fasta sequence of both species to be compared. 1 file per species
input_Folder=$(pwd)

orthofinder -f $input_Folder -og # -og: Stop after inferring orthogroups

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################