######################################################
#### This script needs to run once per species   #####
######################################################

######################### Script will do ########################################## 
# * STEP 1: wallace                                                               #
###################################################################################

### Needed files
# 1- Target species protein fasta
# 2- Outgroup species protein fasta

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

# simplify header 
sed -i'' 's/ uncharacterized.*//' "$Abb_is"_protein.fasta
sed -i'' 's/ LOW.*//' "$Abb_is"_protein.fasta

echo $(pwd)

# orthofinder do not complete on wallace. It stay running for a while, so best to do in locally in my mac.

######################### Script will do ########################################## 
# * STEP 2: mac                                                                   #
###################################################################################

# copy both fasta to mac
cd /Users/danilo/Dropbox/Backup/MacBookPro/PostDoc/Stukenbrock/Project/2020_postdoc_main/03_results/20spp_PA_MajorGenomeAlignment/Wallace_run1/2020_POSTDOC_MAIN/Podospora/0_data/5_orthofinder

Species_is="Verticillium"

rsync -v -a pereira@wallace:"/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/8_select_orthologs/*fasta" .

conda activate env_orthofinder

# run it. after -f it should be a folder containing the protein fasta sequence of both species to be compared. 1 file per species
input_Folder=$(pwd)

orthofinder -f $input_Folder -og

rsync -v -a OrthoFinder pereira@wallace:"/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/8_select_orthologs/"

# on wallace do (for script 6)
cp /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/4_BRAKER2/braker/augustus.hints.codingseq $Output_folder/augustus.hints.codingseq.fasta
