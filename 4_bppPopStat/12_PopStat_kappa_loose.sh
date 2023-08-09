################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
###############################################################################
#   This script will take the each gene alignment and calculate kappa         #
#                                                                             #
# Files needed:                                                               #
#             * Output from script 10_Gap_threshold.v2.sh                     #
#             * parameter file to get kappa (provided "EstimateKAPPA.bpp")    #
#                                                                             #
###############################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

# manual http://biopp.univ-montp2.fr/manual/html/bppsuite/2.4.0/bppsuite.html#Available-statistics
# go to 4.11 BppPopStats: Bio++ Population Genetics Statistics

# Set variables:
Species_is="Verticillium"
Abb_is="Vdahliae"

# folders and files
input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/4_final_gene_set"
output_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/1_loose_kappa"
param_FILE="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/1_EstimateKAPPA.bpp"

# wallace
echo -e '#!/bin/bash' > script.popstat.sh
IFS="|"
while read line
do

linearray=( $line )
outgroup_ID=${linearray[1]}
gene_ID=${linearray[0]}
#echo $outgroup_ID

# param refers to a parameter file needed by bpppopstat (an example param file for this step is provided 'EstimateKAPPA.bpp')
echo "/home/pereira/software/BPP_danilo/bin/bpppopstats param=$param_FILE isolate_fasta=NT_$gene_ID outgroup_in_fasta=$outgroup_ID input_FOLDER=$input_FOLDER output_FOLDER=$output_FOLDER" >> script.popstat.sh

done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/popstat_gene_outgroup.txt



sbatch --job-name=VD_scp12 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=FAIL --mail-user=pereira@evolbio.mpg.de --partition=standard < script.popstat.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
