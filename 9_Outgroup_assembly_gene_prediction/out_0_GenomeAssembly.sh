################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
######################################################################################
#   This script will take raw illumina reads and perform de novo assembly            #
#                                                                                    #
# Files needed:                                                                      #
#             * SRR code of illumina reads                                           #
#                                                                                    #
######################################################################################
#                                                                       
#
#
######################################################################################
Species_is="Verticillium"

mkdir /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/
mkdir /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/0_rawfastq
mkdir /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/1_trimmed
mkdir /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/2_spades
mkdir /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/3_quast

# create script
# from download to assembly
TRIMM="/home/pereira/software/Trimmomatic-0.39"
IlluminaAdaptor="/home/pereira/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
TrimmedFolder="/home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/1_trimmed"

echo -e '#!/bin/bash' > spades.sh
echo "/home/pereira/software/sratoolkit.2.11.0-centos_linux64/bin/prefetch SRR2012746" >> spades.sh
echo "/home/pereira/software/sratoolkit.2.11.0-centos_linux64/bin/vdb-validate SRR2012746/SRR2012746.sra" >> spades.sh
echo "/home/pereira/software/sratoolkit.2.11.0-centos_linux64/bin/fasterq-dump -S SRR2012746/SRR2012746.sra -O /home/pereira/2020_POSTDOC_MAIN/outgroups/Verticillium/0_rawfastq -t /home/pereira/2020_POSTDOC_MAIN/outgroups/Verticillium/0_rawfastq/temp -e 4" >> spades.sh
echo "gzip SRR2012746_1.fastq" >> spades.sh
echo "gzip SRR2012746_2.fastq" >> spades.sh
echo "java -jar $TRIMM/trimmomatic-0.39.jar PE -threads 4 /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/0_rawfastq/SRR2012746_1.fastq.gz /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/0_rawfastq/SRR2012746_2.fastq.gz $TrimmedFolder/fp_VD_outgroup.fastq.gz $TrimmedFolder/fu_VD_outgroup.fastq.gz $TrimmedFolder/rp_VD_outgroup.fastq.gz $TrimmedFolder/ru_VD_outgroup.fastq.gz ILLUMINACLIP:$IlluminaAdaptor:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:10 MINLEN:50" >> spades.sh
echo "/home/pereira/software/SPAdes-3.14.1-Linux/bin/spades.py -t 4 --careful -k 21,33,55,65,69,87,93 -o /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/2_spades --pe1-1 $TrimmedFolder/fp_VD_outgroup.fastq.gz --pe1-2 $TrimmedFolder/rp_VD_outgroup.fastq.gz" >> spades.sh
echo "/data/biosoftware/quast/quast-4.6.3/quast.py /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/2_spades/scaffolds.fasta -t 4 -o /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/3_quast" >> spades.sh

sbatch --job-name=VD_spa --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=4 --time=96:00:00 --mem=120G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=standard spades.sh


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################


