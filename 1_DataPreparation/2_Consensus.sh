################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#####################################################################
# Script used to extract a genome fasta per sample in the vcf file  #
#                                                                   #
# Files needed:                                                     #
#               * Hard-filtered VCF file                            #
#                                                                   #
#####################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################
#
### need change per species
Species_is="Verticillium"
REF="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/0_references/0_genome_DNA/Vdahliae_genomic"

# common variables
BASE="/home/pereira/2020_POSTDOC_MAIN"

VCF_folder="$BASE/$Species_is/1_genome_alignment/4_genotyped_vcf"
VCF_hardFiltered=$( find $VCF_folder -name "*HardFilter.snpONLY.vcf.gz" ) # will get the HardFiltered VCF with path
VCF_BI_nomissing=$( echo "$VCF_hardFiltered" | sed -e "s/.vcf.gz//g" )

Consensus_output="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/3_consensus"

# copy reference genome and index file
cp "$REF"* /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/2_reference_files

cp "$VCF_BI_nomissing.NOmissing.BIallelic.vcf.gz" /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/1_working_files/

# go to script submission folder
cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/3_consensus

# keep SNPs without missing information and BI allelic SNPs. This is necessary or bcftools consensus will not work
echo -e '#!/bin/bash' > 2_script.sh
echo "vcftools --gzvcf $VCF_hardFiltered --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > $VCF_BI_nomissing.NOmissing.BIallelic.vcf.gz" >> 2_script.sh

# index file is necessary for bcftools consensus
echo "tabix -p vcf $VCF_BI_nomissing.NOmissing.BIallelic.vcf.gz" >> 2_script.sh

# extract for each isolate a complete genome fasta file 

# in this loop for mac, it might be needed to use bash `exec bash`. Then to go back `exec zsh` (This won't affect new terminal windows or anything, but it's convenient.)
IFS="|"
while read line;
do
    #echo $line
    linearray=( $line )
    Sample_Name=${linearray[9]}
    echo $Sample_Name
    echo "bcftools consensus --fasta-ref "$REF".fasta "$VCF_BI_nomissing".NOmissing.BIallelic.vcf.gz --sample $Sample_Name -o $Consensus_output/$Sample_Name.genome.fasta" >> 2_script.sh
done < $BASE/$Species_is/$Species_is.global.ratesOfadaptation.txt


sbatch --job-name=VD_spt2 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=24:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=fast 2_script.sh



######################################################################################################################
#                                                      END                                                           #
######################################################################################################################