# manual http://biopp.univ-montp2.fr/manual/html/bppsuite/2.4.0/bppsuite.html#Available-statistics
# go to 4.11 BppPopStats: Bio++ Population Genetics Statistics

# dN_dS
# For codon sequences only. Build the consensus sequence of both ingroup and outgroup alignments and fit a Yang and Nielsen model of 
# codon sequence evolution with a maximum likelihood approach. Reports the estimated parameters ***omega*** (dN / dS ratio) and 
# kappa (transitions / transversions ratio), as well as the divergence between the two sequences.

# Installing this shit is hell!!!

# run examples
# {program} parameter1=value1 parameter2=value2 ... parameterN=valueN
# or
# {program} param=option_file

# run 5 - After Julien feedback on March 30

# Necessary variables:
# 1 - $(input_FOLDER)
# 2 - $(output_FOLDER)
# 3 - $(param_FILE)
# 4 - $targetspp_ID
# 4 - $outgroup_ID

Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/1_loose_kappa

# folders and files
cp /home/pereira/2020_POSTDOC_MAIN/Penicillium/3_analysis/11_popstat/0_param/1_EstimateKAPPA.bpp /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/

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

# Wallace
echo "/home/pereira/software/BPP_danilo/bin/bpppopstats param=$param_FILE isolate_fasta=NT_$gene_ID outgroup_in_fasta=$outgroup_ID input_FOLDER=$input_FOLDER output_FOLDER=$output_FOLDER" >> script.popstat.sh

done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/popstat_gene_outgroup.txt



sbatch --job-name=VD_scp12 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=FAIL --mail-user=pereira@evolbio.mpg.de --partition=standard < script.popstat.sh

