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

# per species
Species_is="Verticillium"
Abb_is="Vdahliae"

cd /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/2_fixed_kappa

# kappa is
cat /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/0_scripts/13_popstat/1_loose_kappa/Vdahliae_mediankappa.txt

# change kappa in the bpp file
nano /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/2_fixedKAPPA.bpp


# folders and files
input_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/10_alignment/4_final_gene_set"
output_FOLDER="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/3_fixed_kappa"
param_FILE="/home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/2_fixedKAPPA.bpp"

# wallace
COUNTER=0
IFS="|"
while read line
do

linearray=( $line )
outgroup_ID=${linearray[1]}
gene_ID=${linearray[0]}
#echo $outgroup_ID

COUNTER=$((COUNTER + 1)) # creates a counter with increment of +1 every time the loop comes back here
printf -v COUNTER_four "%04d" $COUNTER # make the COUNTER 4-digit length and save to another variable COUNTER_four
echo $COUNTER_four
echo -e '#!/bin/bash' > $COUNTER_four.script_$gene_ID.sh

# Wallace
echo "/home/pereira/software/BPP_danilo/bin/bpppopstats param=$param_FILE isolate_fasta=NT_$gene_ID outgroup_in_fasta=$outgroup_ID input_FOLDER=$input_FOLDER output_FOLDER=$output_FOLDER" >> $COUNTER_four.script_$gene_ID.sh

done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/15_popstat/0_param/popstat_gene_outgroup_final.txt

# tests in my mac took > 7min for each run, so go for arrays

ls *.sh > list_of_jobs.txt
cat list_of_jobs.txt | wc -l

# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e 'module load java/x64/8u121' >> array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_four "%04d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 4-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list_of_jobs.txt | grep "$SLURM_ARRAY_TASK_ID_four.script_*")' >> array_batch.sh
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh

# submit
sbatch --job-name=VD_sct13 --array=1-3958%100 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=18G --error=job_%A_%a.err --output=job_%A_%a.out --mail-type=NONE --mail-user=pereira@evolbio.mpg.de --partition=standard array_batch.sh

