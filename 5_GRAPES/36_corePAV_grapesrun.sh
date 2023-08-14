# 16 JUN 2022
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

# mkdir
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/0_bin
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/1_input
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/2_output
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/0_bin
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/1_input
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/2_output
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/0_bin
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/1_input
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/2_output

# move files to rsync -v -a *.tar.gz sunbo481@nesh-fe.rz.uni-kiel.de:"/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/"

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/0_FullDataset_Bootstrap.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/1_input/ 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/1_input/ && tar -xzf 0_FullDataset_Bootstrap.tar.gz

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/1_FullDataset_Permutation.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/1_input/ 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/1_input/ && tar -xzf 1_FullDataset_Permutation.tar.gz

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/2_FullDataset_indv_cat.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/1_input 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/1_input/ && tar -xzf 2_FullDataset_indv_cat.tar.gz 


######################################
#       Dataset complete             #
######################################

######################################
#            bootstrap               #
######################################
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/0_bin

# folders
input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/1_input/0_FullDataset_Bootstrap"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/1_bootstrap/2_output"


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

echo "sbatch --job-name="$Abb_slurm"_BoAll --array=1-200 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_boot_allgenes.sh" > submit_Boall.sh

bash submit_Boall.sh

######################################
#            permutation             #
######################################
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/1_input/1_FullDataset_Permutation"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_perm_allgenes.txt
sed -i 's/.sfs//' SFS_perm_allgenes.txt
mkdir log

# create array batch file
touch jobs_SFS_perm_allgenes.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model $Best_grapes_model_allgenes -no_div_param" >> jobs_SFS_perm_allgenes.sh
done < SFS_perm_allgenes.txt

cat SFS_perm_allgenes.txt | wc -l # this is the number of jobs in the array 2000

# Create the submission script
echo -e '#!/bin/bash' > array_perm_allgenes.sh
echo -e '# get sample' >> array_perm_allgenes.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_perm_allgenes.sh)' >> array_perm_allgenes.sh
echo -e '# run the command' >> array_perm_allgenes.sh
echo -e '$sample_ID' >> array_perm_allgenes.sh

echo "sbatch --job-name="$Abb_slurm"_PeAll --array=1-2000 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_perm_allgenes.sh" > submit_PeAll.sh

bash submit_PeAll.sh


#######################################################
#            all genes run for the mean               #
#######################################################

cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/1_input/"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/3_singleRun/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_run_for_mean.txt
sed -i 's/.sfs//' SFS_run_for_mean.txt
mkdir log

# create array batch file
touch jobs_SFS_run_for_mean.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model $Best_grapes_model_allgenes -no_div_param" >> jobs_SFS_run_for_mean.sh
done < SFS_run_for_mean.txt

cat SFS_run_for_mean.txt | wc -l # this is the number of jobs in the array 300

# Create the submission script
echo -e '#!/bin/bash' > array_run_for_mean.sh
echo -e '# get sample' >> array_run_for_mean.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_run_for_mean.sh)' >> array_run_for_mean.sh
echo -e '# run the command' >> array_run_for_mean.sh
echo -e '$sample_ID' >> array_run_for_mean.sh

echo "sbatch --job-name="$Abb_slurm"_1st --array=1-2 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_run_for_mean.sh" > submit_1st.sh

bash submit_1st.sh


###########################
# run failed permutations #
###########################
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/8_PAV_vs_Core/7_grapes_run/2_permutation/0_bin

module load nano
nano failed_runs.txt

# extract failed runs from job file
while read line
do
    grep "$line" jobs_SFS_perm_allgenes.sh >> jobsFAILED_SFS_perm_allgenes.sh
done < failed_runs.txt

sed -i 's/no_div_param/no_div_param -nb_rand_start 5/g' jobsFAILED_SFS_perm_allgenes.sh

Number_of_array=$( cat failed_runs.txt | wc -l )

# create submission files
# Create the submission script
echo -e '#!/bin/bash' > array_perm_allgenesFAILED.sh
echo -e '# get sample' >> array_perm_allgenesFAILED.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobsFAILED_SFS_perm_allgenes.sh)' >> array_perm_allgenesFAILED.sh
echo -e '# run the command' >> array_perm_allgenesFAILED.sh
echo -e '$sample_ID' >> array_perm_allgenesFAILED.sh

echo "sbatch --job-name="$Abb_slurm"_perm --array=1-$Number_of_array --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --qos=long --time=96:00:00 --mem=42G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster < array_perm_allgenesFAILED.sh" > submit_PeAll_FAILED.sh

bash submit_PeAll_FAILED.sh
