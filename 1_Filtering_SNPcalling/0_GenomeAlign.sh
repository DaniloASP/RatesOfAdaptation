################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
######################################################################################################################
# Script used to clean and perform sanity checks on the reference genome and annotation file of a reference species  #
#                                                                                                                    #
# Files needed:                                                                                                      #
#               * Reference genome in fasta format                                                                   #
#               * Fasta with protein sequences                                                                       #
#               * GFF from reference species                                                                         #
#               * Raw fastq illumina reads                                                                           #
######################################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################
#
#
#
##################
#  Sanity check  #
##################
#
# If necessary simplify contig names
sed 's/ .*//' GCF_002742065.1_CB0940_V2_genomic.fna > CB0940_V2_genomic.fa # it will remove all characters after a first space

# If necessary check for duplicate entries or pseudo genes in the gff file
agat_sp_keep_longest_isoform.pl -gff Vdahliae_genomic.gff -o Vdahliae_noDup_genomic.gff

/Users/danilo/Documents/Git_Repository/gffread/gffread --no-pseudo Vdahliae_noDup_genomic.gff > Vdahliae_noDup_genomic.checked.gff

# If necessary match the name of genes in the GFF and the name of proteins in the fasta file
grep 'XP_' Vdahliae_genomic.gff | cut -d';' -f1,2 | cut -d'=' -f2,3 | sed 's/cds-//' | sed 's/Parent=rna-//' | sed 's/;/|/' | sed 's/Parent=rna-//' | sort | uniq > Script10_fix.txt

echo -e '#!/bin/bash' > script_sed_name.sh
IFS=$'|'
while read line
    do
    linearray=( $line )
    new_name=${linearray[0]}
    old_name=${linearray[1]}
        #echo $length_id
    #echo ">Ztritici:$CHR_id:1:+:$length_id"
    echo "sed -i '' ""'s/[[:<:]]$old_name[[:>:]]/$new_name/g'"" Vdahliae_noDup_genomic.checked.gff" >> script_sed_name.sh
done < Script10_fix.txt

bash script_sed_name.sh

###########################################
# Trimming and mapping reads to reference #     
###########################################
#
# Before starting the loop do the following:
#       1 - Genome reference fasta indexed (bt2) [bowtie2-build Vdahliae_genomic.fasta Vdahliae_genomic]
#       2 - Genome reference fasta index (.fai) [samtools faidx Vdahliae_genomic.fasta]
#       3 - Genome reference fasta dictionary (dict) [gatk CreateSequenceDictionary -R Vdahliae_genomic.fasta]

# These analysis were performed on Wallace-HPC from the MPI Ploen - Slurm managing system
module load java
module load python/3.7.1

# Set path to softwares and dependencies
TRIMM="/home/pereira/software/Trimmomatic-0.39"
IlluminaAdaptor="/home/pereira/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"

# Set pipeline - Note that I use a global file containing information for each isolate in each row. Different collumns represent vary information as isolate ID, raw reads name etc. This information is standardized across species, thus the arrays will not change.
Species_is="Verticillium"
IFS="|"
while read line
do
n=4
#echo $line
linearray=( $line )
givenRawR1=${linearray[5]}
givenRawR2=${linearray[6]}
givenSampleR1=${linearray[7]}
givenSampleR2=${linearray[8]}
Sample=${linearray[9]}

RawReadsFolder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/1_raw_reads/0_DNA"
rawR1=$( find $RawReadsFolder -name "*${givenRawR1}" )
rawR2=$( find $RawReadsFolder -name "*${givenRawR2}" )

#echo $rawR1

#echo $Sample $rawR1 $rawR2

TrimReadsFolder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/2_trimmed_reads/0_DNA"
trimR1="$TrimReadsFolder/fp_$givenSampleR1"
trimR2="$TrimReadsFolder/rp_$givenSampleR2"
unpR1="$TrimReadsFolder/fu_$givenSampleR1"
unpR2="$TrimReadsFolder/ru_$givenSampleR2"

#echo $trimR1 $trimR2

REF="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/0_references/0_genome_DNA/Vdahliae_genomic"
SamFolder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/1_genome_alignment/0_sam"
BamFolder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/1_genome_alignment/1_bam"
gVCFFolder="/home/pereira/2020_POSTDOC_MAIN/$Species_is/1_genome_alignment/2_gvcf"

echo $Sample

echo -e '#!/bin/bash' > $Sample.sh

echo "java -jar $TRIMM/trimmomatic-0.39.jar PE -threads $n $rawR1 $rawR2 $trimR1 $unpR1 $trimR2 $unpR2 ILLUMINACLIP:$IlluminaAdaptor:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:10 MINLEN:50" >> $Sample.sh

echo "bowtie2 -p $n --very-sensitive-local --rg-id $Sample --rg SM:$Sample -x $REF -1 $trimR1 -2 $trimR2 -S $SamFolder/$Sample.sam" >> $Sample.sh

echo "samtools view -bS $SamFolder/$Sample.sam | samtools sort -o $BamFolder/$Sample.bam" >> $Sample.sh

echo "samtools index $BamFolder/$Sample.bam" >> $Sample.sh

echo "gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$REF".fasta -ploidy 1 --emit-ref-confidence GVCF -I $BamFolder/$Sample.bam -O $gVCFFolder/$Sample.g.vcf" >> $Sample.sh

echo "sbatch --job-name=GATK_VD --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=$n --time=72:00:00 --mem=56G --error=job.%J.err --output=job.%J.out --mail-type=FAIL,END --mail-user=pereira@evolbio.mpg.de --partition=standard $Sample.sh" >> sbatch_GATK

echo "sleep 0.2" >> sbatch_GATK

done < /home/pereira/2020_POSTDOC_MAIN/$Species_is/$Species_is.global.ratesOfadaptation.txt


# submit to the queueing system
sh sbatch_GATK

#########################################
# Combine the previously generated gvcf #     
#########################################
Species_is="Verticillium"

# CombineGVCFs
REF="/home/pereira/2020_POSTDOC_MAIN/$Species_is/0_data/0_references/0_genome_DNA/Vdahliae_genomic"
BASE="/home/pereira/2020_POSTDOC_MAIN"

# The first step is to make a list of all the g.vcf files generated with the HaplotypeCaller.
find $BASE/$Species_is/1_genome_alignment/2_gvcf -name "*.g.vcf" > VD_105.g.vcf.list

# Then, it is possible to join all the files with the tool CombineGVCFs
module load java

echo -e '#!/bin/bash' > run_VD_105.sh
echo "gatk --java-options '-XX:+UseSerialGC' CombineGVCFs -R "$REF".fasta \
-V $BASE/$Species_is/2_scripts/2_SNP_calling/VD_105.g.vcf.list  \
-O $BASE/$Species_is/1_genome_alignment/3_combined_gvcf/VD_105.combined.g.vcf" >> run_VD_105.sh

#####################################################
# run GATK4 GenotypeGVCFs on each file individually #
#####################################################
#
# run GenotypeGVCF
# --allow-old-rms-mapping-quality-annotation-data to allow usage of g.vcf files made with gatk versions older than 4
#echo -e '#!/bin/bash' > run_VD_105.sh
echo "gatk --java-options '-XX:+UseSerialGC' GenotypeGVCFs --allow-old-rms-mapping-quality-annotation-data \
	   -R "$REF".fasta \
	   --verbosity ERROR \
	   -V $BASE/$Species_is/1_genome_alignment/3_combined_gvcf/VD_105.combined.g.vcf \
	   --max-alternate-alleles 2 \
	   -O $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotypedGVC.vcf" >> run_VD_105.sh

###################################
# SelectVariants. INDELs and SNPs #
###################################
echo "gatk SelectVariants \
	    -R "$REF".fasta \
	    -V $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotypedGVC.vcf \
	    --select-type-to-include INDEL \
		--select-type-to-include SNP \
	    -O $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotyped.SNP.and.INDEL.vcf" >> run_VD_105.sh

##################################################
# Set filtering thresholds for future hardfilter #
##################################################
#
AN=5 # number of samples
QUAL=1000.0 # How confident we are that there is some kind of variation at a given site. The variation may be present in one or more samples
QD=5.0 # This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage
MQ=20.0 # It is meant to include the standard deviation of the mapping qualities, allowing to check for the variation in the dataset
ReadPosRankSum_lower=-2.0 # This is the u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors. A value close to zero is best because it indicates there is little difference between the positions of the reference and alternate alleles in the reads.
ReadPosRankSum_upper=2.0
MQRankSum_lower=-2.0 # This is the u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping qualities of the reads supporting the reference allele and the alternate allele. A positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele; a negative value indicates the mapping qualities of the reference allele are higher than those supporting the alternate allele. A value close to zero is best and indicates little difference between the mapping qualities.
MQRankSum_upper=2.0
BaseQRankSum_lower=-2.0 # This variant-level annotation tests compares the base qualities of the data supporting the reference allele with those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference
BaseQRankSum_upper=2.0

echo "gatk --java-options '-XX:+UseSerialGC' VariantFiltration -R "$REF".fasta \
	   	-V $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotyped.SNP.and.INDEL.vcf \
	   	--filter 'AN < $AN' --filter-name 'nSamples' \
	   	--filter 'QD < $QD' --filter-name 'QDFilter' \
	   	--filter 'QUAL < $QUAL' --filter-name 'QualFilter' \
	   	--filter 'MQ < $MQ' --filter-name 'MQ' \
	   	--filter 'ReadPosRankSum < $ReadPosRankSum_lower' --filter-name 'ReadPosRankSum' \
	   	--filter 'ReadPosRankSum > $ReadPosRankSum_upper' --filter-name 'ReadPosRankSum' \
	   	--filter 'MQRankSum < $MQRankSum_lower' --filter-name 'MQRankSum' \
	   	--filter 'MQRankSum > $MQRankSum_upper' --filter-name 'MQRankSum' \
	   	--filter 'BaseQRankSum < $BaseQRankSum_lower' --filter-name 'BaseQRankSum' \
	   	--filter 'BaseQRankSum > $BaseQRankSum_upper' --filter-name 'BaseQRankSum' \
	   	-O $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotyped.SNP.and.INDEL.AN$AN.vcf" >> run_VD_105.sh


###########################################################################################
# hard filtering - remove SNPs that have failled one of the previous quality requirements #
###########################################################################################
#
echo "vcftools --vcf $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotyped.SNP.and.INDEL.AN$AN.vcf --remove-filtered-all --remove-indels --recode --stdout | bgzip -c > $BASE/$Species_is/1_genome_alignment/4_genotyped_vcf/VD_105.genotyped.SNP.and.INDEL.AN$AN.HardFilter.snpONLY.vcf.gz" >> run_VD_105.sh

sbatch --job-name=VD_105 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=72:00:00 --mem=32G --error=job.%J.err --output=job.%J.out --mail-type=END,FAIL --mail-user=pereira@evolbio.mpg.de --partition=standard run_VD_105.sh



######################################################################################################################
#                                                      END                                                           #
######################################################################################################################