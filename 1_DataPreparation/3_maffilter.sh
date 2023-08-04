################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#############################################################################################################################################################################
# Script used to extract a CDS features from all sampled genomes                                                                                                            #
#                                                                                                                                                                           #
# Files needed:                                                                                                                                                             #
# Files required:   * Reference genome, with core chr only, in fasta format and fai [samtools faidx Sn15.2021_genome.CORE.fa] (used by bowtie2, but CORE only)              #
#                   * param file for maffilter                                                                                                                              #
#                   * folder with new genomes                                                                                                                               #
#                   * GFF file. Matching contigs/chromosome with the genome sample. NOTE: Very important to have a standard naming for the CDSs and genes.                  #
#                                                                                                                                                                           #
#############################################################################################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################
#
#######################################################################
#   Match coding needed for maffilter for each codon/exon             #
#######################################################################
# Get chr name and length from CORE genome fasta file (reference fasta)
samtools faidx $CORE_reference_fasta

# In the output fai, first column is the CHR id and second the lenght. Change *.fai file delimter to correct one.
Species_is="Verticillium"
Abb_is="Vdahliae"
REF="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/0_references/0_genome_DNA/Vdahliae_genomic"
GFF_name="Vdahliae_noDup_genomic.checked.gff"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/2_reference_files
cut -d$'\t' -f1,2 $REF.fasta.fai > CHR_length.txt

# As JUN 29 2021, a work around to replacing partial match is to use delimter. GNU linux uses \b, but here in mac \b do not work, instad, it is [[:<:]] and [[:>:]]. So I added these to the chr name to avoid replacement of wrong stuff, but [[:<:]] do not work with >, as [[:<:]]>chr..., so the extra code sed -i '' 's/>>/>/g' *.genome.fasta to fix the issue, for more check https://unix.stackexchange.com/questions/190334/sed-word-boundaries-on-macos

# this loop takes about 1h
Output_folder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/3_consensus"

echo -e '#!/bin/bash' > script_sed_chr_name.sh
IFS=$'\t'
while read line
do
linearray=( $line )
CHR_id=${linearray[0]}
length_id=${linearray[1]}
#echo $length_id
#echo ">Ztritici:$CHR_id:1:+:$length_id"
echo "sed -i 's/\b$CHR_id\b/$Abb_is:$CHR_id:1:+:$length_id/g' $Output_folder/*.genome.fasta" >> script_sed_chr_name.sh
done < CHR_length.txt
#echo "sed -i '' 's/>>/>/g' $Output_folder/*.genome.fasta" >> script_sed_chr_name.sh

sbatch --job-name=VD_s3.1 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=2 --time=24:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=fast script_sed_chr_name.sh


################
# organize GFF #
################
#
# I need to standardize the GFF contigs name to the reference genome used to generate the vcf file.

# replace names
IFS='|'
while read line
do
linearray=( $line )
OLD_name=${linearray[0]}
NEW_name=${linearray[1]}
#echo $NEW_name
echo "sed -i 's/\b$OLD_name\b/$NEW_name/g' $GFF_name" >> script_sed.sh
done < key_NCBIcode_contig.txt 

#########################################################################
# Here I use maffilter to extract CDS/exon from each sample genome file #
#########################################################################

# maffilter Wallace
cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/4_maffilter_cds_extr

input_Folder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/3_consensus"
dependencies_Folder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/4_maffilter/0_files_maffilter"
output_Folder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/4_maffilter/1_maffilter_output"
param_file="param.VD.FEB2022"

# loop
IFS="|"
while read line
do
    linearray=( $line )
    Sample_Name=${linearray[9]}
    #echo $Sample_Name
    echo -e '#!/bin/bash' > $Sample_Name.maffilter.sh
    echo "/home/pereira/software/MafFilter/maffilter param=$dependencies_Folder/$param_file Isolate=$Sample_Name input_Folder=$input_Folder output_Folder=$output_Folder dependencies_Folder=$dependencies_Folder" >> $Sample_Name.maffilter.sh
    echo "sbatch --job-name=VD_s3.3 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=24:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=FAIL --mail-user=pereira@evolbio.mpg.de --partition=fast $Sample_Name.maffilter.sh" >> sbatch_MAFFILTER.sh
    echo "sleep 0.2" >> sbatch_MAFFILTER.sh
done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/$Species_is.global.ratesOfadaptation.txt

# execute with bash sbatch_MAFFILTER.sh


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################