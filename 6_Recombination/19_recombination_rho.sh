###########################################################################################################################
# Danilo Pereira - Feb 2022 - Kiel

# Script to develop pipeline for estimation of recombination map across the genome from population genetic data

# paper: https://academic.oup.com/genetics/article/208/3/1209/6066459#supplementary-data
# In their paper, they tested LDhat and LDhelmet, and sticked to LDhat

# tool: https://github.com/auton1/LDhat

# tool manual: https://github.com/auton1/LDhat/blob/master/manual.pdf

# Script version 0.1 

###########################################################################################################################

### mandatory files and naming #####################################################################
# * A .fai file in $BASE/0_data/0_references/0_genome_DNA/*.fai
# * A VCF ending in "*NOmissing.BIallelic.vcf.gz" in $BASE/1_genome_alignment/4_genotyped_vcf 
# * A fasta reference genome ending in "*_genomic.fasta" in $BASE/0_data/0_references/0_genome_DNA
#####################################################################################################

####################################################
### This part needs to be modified per species   ###
####################################################

Species_is="Zymoseptoria_RUN2"
Abb_is="Zymoseptoria"
Abb_slurm="ZT"
#LK_table_raw="lk_n100_t0.001"
LK_table_raw="lk_n100_t0.01"
number_of_chr=100
ploidy=1

####

################################################################################################
### Further modification is not needed as long as folder structure and naming are the same   ###
################################################################################################

############
### Tools ##
############
LDhat_bin="/data/biosoftware/ldhat/bin"
LFtable_raw_folder="/home/pereira/software/LDhat_LK_table"

#############################
#### Stablish directories ###
#############################
BASE="/home/pereira/2020_POSTDOC_MAIN/$Species_is"

##
mkdir $BASE/3_analysis/0_scripts/15_LDhat
mkdir $BASE/3_analysis/0_scripts/15_LDhat/0_first_script

mkdir $BASE/3_analysis/13_Ldhat
mkdir $BASE/3_analysis/13_Ldhat/1_vcf_files
mkdir $BASE/3_analysis/13_Ldhat/2_sample_fasta
mkdir $BASE/3_analysis/13_Ldhat/3_fasta_concatenated
mkdir $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/complete
mkdir $BASE/3_analysis/13_Ldhat/4_convert
mkdir $BASE/3_analysis/13_Ldhat/5_interval
mkdir $BASE/3_analysis/13_Ldhat/6_stat


for repetition in {1..9}
do
mkdir $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_$repetition
mkdir $BASE/3_analysis/13_Ldhat/4_convert/rep_$repetition
mkdir $BASE/3_analysis/13_Ldhat/5_interval/rep_$repetition
mkdir $BASE/3_analysis/13_Ldhat/6_stat/rep_$repetition

done



# go to cd $BASE/3_analysis/0_scripts/15_LDhat/0_first_script


##########################################
### Part 1: File preparation for LDhat ###
##########################################

# It is necessary to extract individual chromosomes or contigs

# get contig names
cp $BASE/0_data/0_references/0_genome_DNA/*.fai $BASE/3_analysis/0_scripts/15_LDhat/0_first_script
cd $BASE/3_analysis/0_scripts/15_LDhat/0_first_script
#cut -f1 -d$'\t' *fai > contig_list.txt
#cut -f2 -d$'\t' *fai > contig_length.txt
cut -f1,2 -d$'\t' *fai > contig_list_length.txt


# VCF input file
# if "*NOmissing.BIallelic.vcf.gz" not there, use vcftools
# cd $BASE/1_genome_alignment/4_genotyped_vcf
# input_VCF=$( find $BASE/1_genome_alignment/4_genotyped_vcf -name "*.HardFilter.snpONLY.vcf.gz" )
# echo -e '#!/bin/bash' > 3_vcftools.sh
# echo "vcftools --gzvcf $input_VCF --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > $BASE/1_genome_alignment/4_genotyped_vcf/ZT_485.genotyped.SNP.and.INDEL.AN5.HardFilter.snpONLY.NOmissing.BIallelic.vcf.gz" >> 3_vcftools.sh
# sbatch --job-name=filter --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=12:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=pereira@evolbio.mpg.de --partition=fast 3_vcftools.sh
# cat out.log
# cd $BASE/3_analysis/0_scripts/15_LDhat/0_first_script

input_VCF=$( find $BASE/1_genome_alignment/4_genotyped_vcf -name "*NOmissing.BIallelic.vcf.gz" )

# reference genome
Reference_genome=$( find $BASE/0_data/0_references/0_genome_DNA -name "*.fa" )

# get number of isolates in vcf. This is a species with n>100, so number_of_chr is fixed to 100 because of subsampling
#number_of_chr=$( bcftools query -l $input_VCF | wc -l )

# sanity check
echo "$input_VCF" >> $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/troubleshooting.txt
echo "$Reference_genome" >> $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/troubleshooting.txt
echo "$number_of_chr" >> $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/troubleshooting.txt


###
# loop to extract each contig from vcf, and convert each isolate from each contig into one fasta
###
rm $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/*fasta # remove any concatenated file there

IFS=$'\t'
while read line_contig
    do
    linearray_1=( $line_contig )
    Contig_name=${linearray_1[0]}
    Contig_length=${linearray_1[1]}

    echo -e '#!/bin/bash' > 0_$Contig_name.LD_hat.sh

    #echo $Contig_name
    #echo $Contig_length

    # extract a vcf per contig
    echo "vcftools --gzvcf $input_VCF --min-alleles 2 --max-alleles 2 --chr $Contig_name --max-missing 1 --out $BASE/3_analysis/13_Ldhat/1_vcf_files/$Abb_is.chr.$Contig_name --recode --stdout | bgzip -c > $BASE/3_analysis/13_Ldhat/1_vcf_files/$Abb_is.chr.$Contig_name.noMiss.BI.vcf.gz" >> 0_$Contig_name.LD_hat.sh
    echo "tabix -p vcf $BASE/3_analysis/13_Ldhat/1_vcf_files/$Abb_is.chr.$Contig_name.noMiss.BI.vcf.gz" >> 0_$Contig_name.LD_hat.sh

    # Convert the alignment to fasta
    IFS='|'

    while read line_sample
        do
        linearray_2=( $line_sample )
        Sample_name=${linearray_2[9]}
        #echo $Sample_name
        echo "samtools faidx $Reference_genome $Contig_name | bcftools consensus $BASE/3_analysis/13_Ldhat/1_vcf_files/$Abb_is.chr.$Contig_name.noMiss.BI.vcf.gz --sample $Sample_name -o $BASE/3_analysis/13_Ldhat/2_sample_fasta/$Sample_name.chr.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh

        # rename sequence of contig inside the fasta 
        echo "seqtk rename $BASE/3_analysis/13_Ldhat/2_sample_fasta/$Sample_name.chr.$Contig_name.fasta "$Sample_name.chr.$Contig_name" > $BASE/3_analysis/13_Ldhat/2_sample_fasta/temp_$Sample_name.chr.$Contig_name.fasta && rm $BASE/3_analysis/13_Ldhat/2_sample_fasta/$Sample_name.chr.$Contig_name.fasta && mv $BASE/3_analysis/13_Ldhat/2_sample_fasta/temp_$Sample_name.chr.$Contig_name.fasta $BASE/3_analysis/13_Ldhat/2_sample_fasta/$Sample_name.chr.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh

        # concatenate fasta into single file
        echo "cat $BASE/3_analysis/13_Ldhat/2_sample_fasta/$Sample_name.chr.$Contig_name.fasta >> $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/complete/$Species_is.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh
    done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/$Abb_is.global.ratesOfadaptation.txt

    IFS=$'\t'

    # this part is excluse for species that have more than n>100

    ### rep_1
    #fasta-subsample $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/complete/$Species_is.$Contig_name.fasta $number_of_chr -seed 47 > $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_1/rep_1.$Species_is.$Contig_name.fasta

    # concatenate variables and create the first line needed by LDhat convert
    #First_line="${number_of_chr} ${Contig_length} ${ploidy}"
    #echo $First_line

    # add first line to each concatenated chr fasta file
    #echo "sed -i '1i$First_line' $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_1/rep_1.$Species_is.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh

    for rep_num in {1..9}
    do
        #echo $rep_num
        echo "fasta-subsample $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/complete/$Species_is.$Contig_name.fasta $number_of_chr -seed 3$rep_num > $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh

        # concatenate variables and create the first line needed by LDhat convert
        First_line="${number_of_chr} ${Contig_length} ${ploidy}"
        #echo $First_line

        # add first line to each concatenated chr fasta file
        echo "sed -i '1i$First_line' $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.fasta" >> 0_$Contig_name.LD_hat.sh


    done

    echo "sbatch --job-name="$Abb_slurm"_rec1 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=FAIL,END --mail-user=pereira@evolbio.mpg.de --partition=standard 0_$Contig_name.LD_hat.sh" >> 0_sbatch_LDhat.sh


done < $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/contig_list_length.txt



# check if all is good before submission
cat $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/troubleshooting.txt

# bash 0_sbatch_LDhat.sh

###########################
### Part 2: Using LDhat ###
###########################

# * $LDhat_bin/lkgen: If LK table is not at the proper number of chromosomes it need to be downsampled. ALL SPP NEED LD TABLE DOWNSAMPLE.
# * $LDhat_bin/convert: Will take the concatenated fast contigs previously generated and make the input for next step
# * $LDhat_bin/interval: Will use output from convert + LK table and estimate rho
# * $LDhat_bin/stat: Will summarize the results and the analyse is done

### lkgen. No need to downsample, table is already n=100
LDhat_bin="/data/biosoftware/ldhat/bin"
LFtable_raw_folder="/home/pereira/software/LDhat_LK_table"

### convert && interval && stat
IFS=$'\t'
while read line_contig
    do
    linearray_1=( $line_contig )
    Contig_name=${linearray_1[0]}
    Contig_length=${linearray_1[1]}

    echo -e '#!/bin/bash' >> 1_$Species_is.$Contig_name.sh

    for rep_num in {1..9}
    do

        echo "$LDhat_bin/convert -seq $BASE/3_analysis/13_Ldhat/3_fasta_concatenated/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.fasta -prefix $BASE/3_analysis/13_Ldhat/4_convert/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name." >> 1_$Species_is.$Contig_name.sh

        echo "$LDhat_bin/interval -seq $BASE/3_analysis/13_Ldhat/4_convert/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.sites.txt -loc $BASE/3_analysis/13_Ldhat/4_convert/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.locs.txt -lk $LFtable_raw_folder/$LK_table_raw -its 10000000 -samp 5000 -bpen 20 -prefix $BASE/3_analysis/13_Ldhat/5_interval/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name." >> 1_$Species_is.$Contig_name.sh

        echo "$LDhat_bin/stat -input $BASE/3_analysis/13_Ldhat/5_interval/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.rates.txt -burn 1000 -loc $BASE/3_analysis/13_Ldhat/4_convert/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name.locs.txt -prefix $BASE/3_analysis/13_Ldhat/6_stat/rep_$rep_num/rep_$rep_num.$Species_is.$Contig_name. # recomended burn 100000 but I get an error, trying with 1000 works" >> 1_$Species_is.$Contig_name.sh
    
    done


    echo "sbatch --job-name="$Abb_slurm"_rec1 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=24G --error=job.%J.err --output=job.%J.out --mail-type=FAIL,END --mail-user=pereira@evolbio.mpg.de --partition=standard 1_$Species_is.$Contig_name.sh" >> 1_sbatch_LDhat.sh

    echo "sleep 0.2" >> 1_sbatch_LDhat.sh

done < $BASE/3_analysis/0_scripts/15_LDhat/0_first_script/contig_list_length.txt


# bash 1_sbatch_LDhat.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

