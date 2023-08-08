################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#########################################################################################################################################################################################################
# This script will align the outgroup protein sequence to the ingroup alignemnt (the ingroup comes from VCF thus is aligned), avoiding gaps in the ingroup.                                             #
#                                                                                                                                                                                                       #
# Files needed:                                                                                                                                                                                         #
#             * A folder containing a key file                                                                                                                                                          #                             
#             * A folder containing multiple fasta files for each gene                                                                                                                                  #
#             * A folder containing individual target species fasta genes                                                                                                                               #
#                                                                                                                                                                                                       #
# Note for next step in the pipeline: Gaps are coded as "-", and bpppop cannot read it. Thus, replace by "N" using "sed -i '' 's/-/N/g' *.fasta"                                                        #
#                                                                                                                                                                                                       #    
# IMPORTANT: USE BASH AND NOT ZSH. THIS SCRIPT USES ARRAY (while loop). CHANGE TO BASH USING `exec bash`                                                                                                #
#                                                                                                                                                                                                       #                                   
# IMPORTANT: There are degrees of optmization and speed for MACSE. In this script I use a 'relatively quick optimization' given by -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2    #
#                                                                                                                                                                                                       #
# IMPORTANT: This script uses array on wallace for running > 7000 jobs.                                                                                                                                 #
#                                                                                                                                                                                                       #
#########################################################################################################################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

## change needed per species
Species_is="Verticillium"
Abb_is="Vdahliae"

# go to cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/10_alignment
## fixed
InputFolder_ingroup="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/4_maffilter/6_split_into_genes_fasta.v2"
InputFolder_outgroup="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/9_outgroup_genes"
OutputFolder_mafft="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/1_mafft_output"
OutputFolder_macse="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/2_macse_output"


# Build scripts for alingment 
module load java/x64/8u121
COUNTER=0
IFS=$'\t'
while read -r line
do
    linearray=( $line )
    Outgroup_gene=${linearray[1]}
    Ingroup_gene=${linearray[0]}
    
    # Wallace
    COUNTER=$((COUNTER + 1)) # creates a counter with increment of +1 every time the loop comes back here
    printf -v COUNTER_four "%04d" $COUNTER # make the COUNTER 4-digit length and save to another variable COUNTER_four
    echo $COUNTER_four
    echo -e '#!/bin/bash' > $COUNTER_four.script_$Ingroup_gene.sh

    # MAFFT
    echo "/home/pereira/software/mafft-7.475/mafft --add $InputFolder_outgroup/$Outgroup_gene.fasta --keeplength $InputFolder_ingroup/$Ingroup_gene.fasta > $OutputFolder_mafft/$Ingroup_gene.fasta" >> $COUNTER_four.script_$Ingroup_gene.sh
    
    # MACSE
    echo "java -jar /home/pereira/software/macse/macse_v2.05.jar -prog refineAlignment -align $OutputFolder_mafft/$Ingroup_gene.fasta -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2 -out_NT $OutputFolder_macse/1_NT/NT_$Ingroup_gene.fasta -out_AA $OutputFolder_macse/2_AA/AA_$Ingroup_gene.fasta" >> $COUNTER_four.script_$Ingroup_gene.sh

    # Wallace
    #echo "sbatch --job-name=MAFMACSE --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=BEGIN,END,FAIL --mail-user=pereira@evolbio.mpg.de --partition=standard script_$Ingroup_gene.sh" >> sbatch_MAFFTMACSE.shls

    #echo "sleep 1" >> sbatch_MAFFTMACSE.sh

done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/8_select_orthologs/Species_OUT_orthologs_keyOnly.txt


# Since at this step there are thousands of genes, meaning, thousands of jobs required, it's best to use array to not overload the submission queue.
# create a list of files necessary for the array extract names
ls *.sh > list_of_jobs.txt
cat list_of_jobs.txt | wc -l # add this number to the upper array value for sbatch: 7544 rows

# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e 'module load java/x64/8u121' >> array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_four "%04d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 4-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list_of_jobs.txt | grep "$SLURM_ARRAY_TASK_ID_four.script_*")' >> array_batch.sh
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh



sbatch --job-name=VD_s7 --array=1-7544 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=24G --error=MAFMACSE_%A_%a.err --output=MAFMACSE_%A_%a.out --mail-type=FAIL --mail-user=pereira@evolbio.mpg.de --partition=standard array_batch.sh



######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

