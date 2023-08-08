########################################################################################
#### This script is intended to copy final set of genes, replace gapes as - to N   #####
########################################################################################

Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/1_loose_kappa

### no modification needed
#

input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/1_loose_kappa"
mkdir /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/1.1_copy_log
copy_folder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/1.1_copy_log"

echo -e '#!/bin/bash' > 12_1_stopcodon_filter.sh
echo "input_FOLDER=$input_FOLDER" >> 12_1_stopcodon_filter.sh
echo "Species_is=$Species_is" >> 12_1_stopcodon_filter.sh

echo "cp $input_FOLDER/* $copy_folder/" >> 12_1_stopcodon_filter.sh
#cp $copy_folder/* $input_FOLDER/

# check number of stop codons, if different then 1, remove gene from next run.
# find genes without stopcodon
echo 'grep -L '"'Info:'"' $input_FOLDER/*.log > $input_FOLDER/0_No_stopcodon.txt' >> 12_1_stopcodon_filter.sh

# find genes with multiple stops
echo 'grep '"'Info:'"' $input_FOLDER/*.log > $input_FOLDER/1_multiple_stopcodon_list.txt' >> 12_1_stopcodon_filter.sh
echo 'sed -i '"'/discarded 1 sites/d'"' $input_FOLDER/1_multiple_stopcodon_list.txt' >> 12_1_stopcodon_filter.sh

# merge and clean
echo 'cat $input_FOLDER/0_No_stopcodon.txt $input_FOLDER/1_multiple_stopcodon_list.txt > $input_FOLDER/2_remove_genes.txt' >> 12_1_stopcodon_filter.sh
echo 'cat $input_FOLDER/2_remove_genes.txt | cut -d'"':'"' -f1 > $input_FOLDER/3_remove_genes.clean.txt' >> 12_1_stopcodon_filter.sh

# remove logs from folder
echo 'while read line' >> 12_1_stopcodon_filter.sh
echo 'do' >> 12_1_stopcodon_filter.sh
echo 'rm $line' >> 12_1_stopcodon_filter.sh
echo 'done < $input_FOLDER/3_remove_genes.clean.txt' >> 12_1_stopcodon_filter.sh

echo 'grep '"'Kappa'"' $input_FOLDER/*log | cut -d'"'/'"' -f9 | sed '"'s/.log:/|/'"'g | sed '"'s/ = /|/'"'g > $input_FOLDER/4_kappa_final.txt' >> 12_1_stopcodon_filter.sh

# match names
echo 'cat $input_FOLDER/3_remove_genes.clean.txt | cut -d'"'/'"' -f9 | sed '"'s/NT_//'"' | sed '"'s/.log//'"' > $input_FOLDER/4_remove_multiple_stopcodons_clean.txt' >> 12_1_stopcodon_filter.sh

# start removing genes from the popstat final list of genes
# option 1. Remove from file 2 names in file 1. from: https://stackoverflow.com/questions/57038167/how-to-remove-lines-based-on-another-file
echo 'grep -wvf $input_FOLDER/4_remove_multiple_stopcodons_clean.txt popstat_gene_outgroup.txt > popstat_gene_outgroup_final.txt' >> 12_1_stopcodon_filter.sh

echo 'cp popstat_gene_outgroup_final.txt /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param' >> 12_1_stopcodon_filter.sh

echo 'cp /home/pereira/2020_POSTDOC_MAIN/Neurospora/3_analysis/11_popstat/0_param/2_fixedKAPPA.bpp /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param' >> 12_1_stopcodon_filter.sh

bash 12_1_stopcodon_filter.sh

# kappa for final gene set is in $input_FOLDER/4_kappa_final.txt
