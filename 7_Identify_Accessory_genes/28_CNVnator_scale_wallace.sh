################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
########################################################################################
#   This script will run CNVnator to identify structural variations across the genome  #
#                                                                                      #
# Files needed:                                                                        #
#             * bam file for each isolate                                              #
#                                                                                      #
########################################################################################
#
#
#
########################################################################################
# 
# JUL 2022
# CNVnator

# go to folder for submission
cd /home/pereira/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp


IFS="|"
while read line
do
linearray=( $line )
Species_is=${linearray[1]} # species name and also folder name
Path_to_reference=${linearray[2]} # full path to reference fasta
Reference_genome=${linearray[3]} # name of reference fasta file
ABB_spp=${linearray[4]} # two letter abbreviation for each species
Global_file=${linearray[5]} # name at the beggining of global file

BASE="/home/pereira/2020_POSTDOC_MAIN/$Species_is"

echo $Path_to_reference/$Reference_genome

# make dirs
mkdir $BASE/3_analysis/15_CNVnator
mkdir $BASE/3_analysis/15_CNVnator/0_bin
mkdir $BASE/3_analysis/15_CNVnator/1_input
mkdir $BASE/3_analysis/15_CNVnator/2_output

# set folders
bin_folder="$BASE/3_analysis/15_CNVnator/0_bin"
input_folder="$BASE/1_genome_alignment/1_bam"
output_folder="$BASE/3_analysis/15_CNVnator/2_output"

# go to working folder
cd $bin_folder
cp "$Path_to_reference/$Reference_genome" $bin_folder

# split reference genome into contigs or chr
/home/pereira/software/seqkit split --by-id $Reference_genome
mv *.split 0_chr_folder

# set folder
chr_folder="$BASE/3_analysis/15_CNVnator/0_bin/0_chr_folder"

# get contigs
species_contig=$( grep '>' "$Path_to_reference/$Reference_genome" | sed 's/>//g' | awk '{print}' ORS=' ' )

# Wallace

echo -e '#!/bin/bash' > $ABB_spp.script.sh

    while read line
    do
    linearray=( $line )
    Sample=${linearray[9]}

    echo $Sample

    ## Step 1: EXTRACTING READ MAPPING FROM BAM/SAM FILES
    echo "/data/biosoftware/cnvnator/CNVnator/cnvnator -root $output_folder/$Sample.root -chrom "$species_contig"-tree $input_folder/$Sample.bam" >> $ABB_spp.script.sh

    ## Step 2: GENERATING A READ DEPTH HISTOGRAM
    echo "/data/biosoftware/cnvnator/CNVnator/cnvnator -root $output_folder/$Sample.root -his 100 -d $chr_folder" >> $ABB_spp.script.sh

    ## Step 3: CALCULATING STATISTICS
    echo "/data/biosoftware/cnvnator/CNVnator/cnvnator -root $output_folder/$Sample.root -stat 100" >> $ABB_spp.script.sh

    ## Step 4: RD SIGNAL PARTITIONING
    echo "/data/biosoftware/cnvnator/CNVnator/cnvnator -root $output_folder/$Sample.root -partition 100" >> $ABB_spp.script.sh

    ## Step 5: CNV CALLING
    echo "/data/biosoftware/cnvnator/CNVnator/cnvnator -root $output_folder/$Sample.root -call 100 > $output_folder/$Sample.CNV.output" >> $ABB_spp.script.sh


    done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/$Global_file.global.ratesOfadaptation.txt

# submit  
# create a list of files necessary for the array extract names

echo "cd $BASE/3_analysis/15_CNVnator/0_bin" >> /home/pereira/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp/sbatch_submit_CNV.sh
echo "mkdir log" >> /home/pereira/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp/sbatch_submit_CNV.sh
echo "sbatch --job-name="$ABB_spp"_cnv --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=24:00:00 --mem=24G --error=log/job.%J.err --output=log/job.%J.out --mail-type=FAIL --mail-user=pereira@evolbio.mpg.de --partition=global $ABB_spp.script.sh" >> /home/pereira/2020_POSTDOC_MAIN/0_CommonAnalysis/3_CNVnator_allspp/sbatch_submit_CNV.sh


done < /home/pereira/2020_POSTDOC_MAIN/1_log_file_20spp.txt 


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

