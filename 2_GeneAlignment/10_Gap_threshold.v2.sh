################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################################
# This script is intended to quantify gaps within alignments                       #
#                                                                                  #
# Files needed:                                                                    #
#             * Folder with MACSE output alignment files (after initial filters)   #
#                                                                                  #
#                                                                                  #
####################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/12_filtering

nano 10.1_Gap_treshold.py

module load python/python-anaconda-5.0.1 
source activate Danilo_postdoc

input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/3_filtering"
#input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/12_filtering/1_quantification/test_fasta"

#rm GAP_inside_fastaNAME.txt
#rm GAP_inside_seqkit.txt
#rm GAP_inside_seqNAME.txt
#rm script_GAP_inside.sh
touch GAP_inside_fastaNAME.txt # will contain information on the fasta file name
touch GAP_inside_seqNAME.txt # will contain information on the sequence name
touch GAP_inside_seqkit.txt # will contain statistics on the number of gaps and length of sequence

# If this loop is executed to echo into a file, it will be too big for Slurm queueing system (7mb), so it will be split like this
echo -e '#!/bin/bash' > script_GAP_inside.sh
echo "Species_is=$Species_is" >> script_GAP_inside.sh
echo "input_FOLDER=$input_FOLDER" >> script_GAP_inside.sh

echo 'for file in $input_FOLDER/*.fasta' >> script_GAP_inside.sh
echo 'do' >> script_GAP_inside.sh
echo 'basename $file' >> script_GAP_inside.sh
echo 'basename "$file" >> GAP_inside_fastaNAME.txt' >> script_GAP_inside.sh
echo 'basename "$file" >> GAP_inside_fastaNAME.txt' >> script_GAP_inside.sh
echo 'awk '"'NR==1{print}'"' $file >> GAP_inside_seqNAME.txt' >> script_GAP_inside.sh # get info on line 1
echo 'awk '"'NR==1{print}'"' $file >> GAP_inside_seqNAME.txt' >> script_GAP_inside.sh
echo 'awk '"'NR==1; NR==2{print}'"' $file | /home/pereira/software/seqkit stats -T -a >> GAP_inside_seqkit.txt' >> script_GAP_inside.sh  # get info on line 1 and 2
echo 'basename "$file" >> GAP_inside_fastaNAME.txt' >> script_GAP_inside.sh
echo 'basename "$file" >> GAP_inside_fastaNAME.txt' >> script_GAP_inside.sh
echo 'tail -n 2 $file | grep '"'>'"' >> GAP_inside_seqNAME.txt' >> script_GAP_inside.sh # get info on second last position
echo 'tail -n 2 $file | grep '"'>'"' >> GAP_inside_seqNAME.txt' >> script_GAP_inside.sh
echo 'tail -n 2 $file | /home/pereira/software/seqkit stats -T -a >> GAP_inside_seqkit.txt' >> script_GAP_inside.sh
echo 'done' >> script_GAP_inside.sh

# put dataframe together
echo "paste -d , GAP_inside_fastaNAME.txt GAP_inside_seqNAME.txt GAP_inside_seqkit.txt > GAP_quantification.txt" >> script_GAP_inside.sh

# get GAP_quantification.txt to mac
#rsync pereira@wallace:"/home/pereira/2020_POSTDOC_MAIN/Penicillium/3_analysis/0_scripts/12_filtering/GAP_quantification.txt" /Users/danilo/Dropbox/Backup/MacBookPro/PostDoc/Stukenbrock/Project/2020_postdoc_main/03_results/13spp_PB_MajorGenomeAlignment/Wallace_run1/2020_POSTDOC_MAIN/Penicillium/0_data/12_filtering

echo "sed -i 's/,/|/g' GAP_quantification.txt" >> script_GAP_inside.sh
echo "sed -i 's/>//g' GAP_quantification.txt" >> script_GAP_inside.sh
echo "sed -i 's/\t/|/g' GAP_quantification.txt" >> script_GAP_inside.sh

echo "python 10.1_Gap_treshold.py" >> script_GAP_inside.sh # execute script

echo "cut -d$'\t' -f1 3_GAP_geneNAME_failed.5.UNIQUEnames.txt > 3_GAP_geneNAME_failed.5.UNIQUEnames.column1.txt" >> script_GAP_inside.sh

# remove genes
echo "cat 3_GAP_geneNAME_failed.5.UNIQUEnames.column1.txt | wc -l > 5_number_of_genes_failed.txt" >> script_GAP_inside.sh
echo "input_FOLDER="$input_FOLDER"" >> script_GAP_inside.sh
echo 'while read line; do rm "$input_FOLDER"/"$line"; done < 3_GAP_geneNAME_failed.5.UNIQUEnames.column1.txt' >> script_GAP_inside.sh

sbatch --job-name=VD_sct10 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=24:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=NONE --mail-user=pereira@evolbio.mpg.de --partition=fast script_GAP_inside.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
