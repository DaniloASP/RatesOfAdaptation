################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
######################################################################################
#   This script will take assembled genomes and perform gene predictions             #
#                                                                                    #
# Files needed:                                                                      #
#             * Fasta file for the genome                                            #
#             * Protein database                                                     #
#                                                                                    #
######################################################################################
#
#
#
######################################################################################


##########################################################
# Super fragmented genoms might be a probmem             #
##########################################################
# BRAKER2 prefers scaffolds larger then at least 1kb, so remove everything bellow that
cd /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/4_BRAKER2/0_ref_files
grep '>' scaffolds.fasta | wc -l
/home/pereira/software/seqkit seq -m 1000 scaffolds.fasta > VD_scaffolds.1kb.fasta
grep '>' VD_scaffolds.1kb.fasta | wc -l


##########################################################
# Prepare fungal database with protein of target species #
##########################################################

# get Fungi protein database from orthodb. Needs to be done once. Then, for each outgroup, add the proteins from target into the "clean" database.
ls /home/pereira/software/OrthoDB_database/Fungi/$Species_is/
cat /home/pereira/software/OrthoDB_database/Fungi/proteins.fasta /home/pereira/software/OrthoDB_database/Fungi/$Species_is/Vdahliae_protein.fasta > /home/pereira/software/OrthoDB_database/Fungi/$Species_is/VD.proteinsDB.fasta

#####################################################################
####         RUN BRAKER2 - with flags from Stauber et al          ###
#####################################################################
# braker2 version 2.1.6
# pipeline used BRAKER with proteins of any evolutionary distance
# what files from the output to use? https://github.com/Gaius-Augustus/BRAKER/issues/194

# SETUP
# OBS: I'm using an installation made by Kristian on wallace.
# Follow instructions in /data/biosoftware/braker2/READMEv2.1.6.txt 

# if white spaces on protein fasta, a warning is issued, but it can continue
#sed 's/ .*//' Sn15.2021_protein.fa > Sn15.2021.nowhite_protein.fa

cd /home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/4_BRAKER2

# references
REF_genome="/home/pereira/2020_POSTDOC_MAIN/outgroups/$Species_is/4_BRAKER2/0_ref_files/VD_scaffolds.1kb"
PROTEIN_database="/home/pereira/software/OrthoDB_database/Fungi/$Species_is/VD.proteinsDB"

# exports and module for wallace
#cp /data/biosoftware/braker2/deps/gm_key_64 $HOME/.gm_key
module load perl/5.26.1
source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
conda activate
module load java/x64/8u121
export PATH=/data/biosoftware/braker2/BRAKER-2.1.6/scripts:$PATH

# create script for wallace using BRAKER-2.1.4
n=8
AUGUSTUS_BASE=/data/biosoftware/braker2/deps/Augustus/

export AUGUSTUS_CONFIG_PATH=$AUGUSTUS_BASE/config/
export AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_BASE/scripts
export AUGUSTUS_BIN_PATH=$AUGUSTUS_BASE/bin
export GENEMARK_PATH=/data/biosoftware/braker2/deps/gmes_linux_64/
export BAMTOOLS_PATH=/data/biosoftware/braker2/deps/bamtools/bin/usr/local/bin/
export DIAMOND_PATH=/data/biosoftware/braker2/deps/diamond/
export BLAST_PATH=/data/biosoftware/braker2/deps/ncbi-blast-2.11.0+/bin/
export PROTHINT_PATH=/data/biosoftware/braker2/deps/ProtHint/bin/
export SAMTOOLS_PATH=/data/biosoftware/braker2/deps/samtools/
export CDBTOOLS_PATH=/data/biosoftware/braker2/deps/cdbfasta/


echo -e '#!/bin/bash' > sbatch_BRAKER2.sh
echo "braker.pl --fungus --cores $n \
--species=VD2_outgroup \
--genome=$REF_genome.fasta \
--prot_seq=$PROTEIN_database.fasta \
--gff3 --alternatives-from-evidence=false" >> sbatch_BRAKER2.sh

sbatch --job-name=VD_break --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --time=72:00:00 --mem=120G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=standard < sbatch_BRAKER2.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

