################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#######################################################################################################################
#   This script will take the SFS files and will calculate the rates of adaptation omegaA, omegaNA and alpha          #
#                                                                                                                     #
# Files needed:                                                                                                       #
#             * SFS files                                                                                             #
#                                                                                                                     #
#######################################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################
# CAU_cluster

# manual https://github.com/BioPP/grapes

####################################################################################

Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=105
Best_grapes_model_allgenes="GammaZero"

# cau_cluster
module load cmake/3.18.4 
module load gcc/10.2.0 
module load eigen/3.3.8
module load gsl/2.5

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/0_FullDataset_Bootstrap.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/1_bootstrap/1_input/ 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/1_bootstrap/1_input/ && tar -xzf 0_FullDataset_Bootstrap.tar.gz

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/1_FullDataset_Permutation.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/2_permutation/1_input/ 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/2_permutation/1_input/ && tar -xzf 1_FullDataset_Permutation.tar.gz

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/2_FullDataset_indv_cat.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/1_input 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/1_input/ && tar -xzf 2_FullDataset_indv_cat.tar.gz 


######################################
#       Dataset complete             #
######################################

######################################
#            bootstrap               #
######################################
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/1_bootstrap/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/1_bootstrap/1_input/0_FullDataset_Bootstrap"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/1_bootstrap/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_boot_allgenes.txt
sed -i 's/.sfs//' SFS_boot_allgenes.txt
mkdir log

# create array batch file
touch jobs_SFS_boot_allgenes.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model $Best_grapes_model_allgenes -no_div_param" >> jobs_SFS_boot_allgenes.sh
done < SFS_boot_allgenes.txt

cat SFS_boot_allgenes.txt | wc -l # this is the number of jobs in the array 300

# Create the submission script
echo -e '#!/bin/bash' > array_boot_allgenes.sh
echo -e '# get sample' >> array_boot_allgenes.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_boot_allgenes.sh)' >> array_boot_allgenes.sh
echo -e '# run the command' >> array_boot_allgenes.sh
echo -e '$sample_ID' >> array_boot_allgenes.sh

echo "sbatch --job-name="$Abb_slurm"_BoAll --array=1-300 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_boot_allgenes.sh" > submit_Boall.sh

bash submit_Boall.sh

######################################
#            permutation             #
######################################


cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/2_permutation/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/2_permutation/1_input/1_FullDataset_Permutation"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/2_permutation/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_perm_allgenes.txt
sed -i 's/.sfs//' SFS_perm_allgenes.txt
mkdir log

# create array batch file
touch jobs_SFS_perm_allgenes.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model $Best_grapes_model_allgenes -no_div_param" >> jobs_SFS_perm_allgenes.sh
done < SFS_perm_allgenes.txt

cat SFS_perm_allgenes.txt | wc -l # this is the number of jobs in the array 3000

# Create the submission script
echo -e '#!/bin/bash' > array_perm_allgenes.sh
echo -e '# get sample' >> array_perm_allgenes.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_perm_allgenes.sh)' >> array_perm_allgenes.sh
echo -e '# run the command' >> array_perm_allgenes.sh
echo -e '$sample_ID' >> array_perm_allgenes.sh

echo "sbatch --job-name="$Abb_slurm"_PeAll --array=1-3000 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_perm_allgenes.sh" > submit_PeAll.sh

bash submit_PeAll.sh


#######################################################
#            all genes run for the mean               #
#######################################################

cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/1_input/"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_run_for_mean.txt
sed -i 's/.sfs//' SFS_run_for_mean.txt
mkdir log

# create array batch file
touch jobs_SFS_run_for_mean.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model all -no_div_param" >> jobs_SFS_run_for_mean.sh
done < SFS_run_for_mean.txt

cat SFS_run_for_mean.txt | wc -l # this is the number of jobs in the array 300

# Create the submission script
echo -e '#!/bin/bash' > array_run_for_mean.sh
echo -e '# get sample' >> array_run_for_mean.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_run_for_mean.sh)' >> array_run_for_mean.sh
echo -e '# run the command' >> array_run_for_mean.sh
echo -e '$sample_ID' >> array_run_for_mean.sh

echo "sbatch --job-name="$Abb_slurm"_1st --array=1-3 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_run_for_mean.sh" > submit_1st.sh

bash submit_1st.sh

###########################################################################
#            all genes run for the SPECIES - get BEST MODEL               #
###########################################################################

cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/1_input/"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/2_output"


echo -e '#!/bin/bash' > all_genes_species.sh
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/1_input/AllCoreGenes.genes.sfs -out /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/7_Gene_function3_core/0_all_genes/3_singleRun/2_output/AllCoreGenes.genes.res -model all -no_div_param" >> all_genes_species.sh

echo "sbatch --job-name="$Abb_slurm"_nonS --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --qos=long --time=96:00:00 --mem=42G --error=log/log.%A.err --output=log/log.%A.out --partition=cluster < all_genes_species.sh" > submit_all_genes_job.sh

bash submit_all_genes_job.sh


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
