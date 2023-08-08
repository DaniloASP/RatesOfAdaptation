################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################################
# This script will annotate functions to the various protein sequences             #
#                                                                                  #
# Files needed:                                                                    #
#             * Target species protein fasta                                       #
#                                                                                  #
# Tools manuals: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html  #
#                                                                                  #
####################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

################
# interproscan #
################
Species_is="Verticillium"
Abb_is="Vdahliae"


# copy data
ls /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan/
cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan

cp /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/8_select_orthologs/*_protein.fasta /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan

# If protein sequences end in "*", it should removed
#sed -i 's/*//g' /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan/Asp_NRRL3357_protein.simpleHeader.fasta

# run.
module load java/x64/11u1
module load perl/5.30.1
module load python/3.7.1
echo -e '#!/bin/bash' > interpro.run.sh
echo "/data/biosoftware/InterProScan/interproscan-5.48-83.0/interproscan.sh -appl SignalP_EUK-4.1 -appl TMHMM-2.0c -dp --goterms --pathways -iprlookup -f TSV, GFF3 -i Vdahliae_protein.fasta" >> interpro.run.sh

# wallace
sbatch --job-name=VD_sc8 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=64G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=standard < interpro.run.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

