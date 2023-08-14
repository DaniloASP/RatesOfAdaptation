# 16 JUN 2022
# CAU_cluster

# manual https://github.com/BioPP/grapes

####################################################################################

Species_is="Verticillium"
Abb_is="Vdahliae"
Abb_slurm="VD"
Num_chr=105
Best_grapes_model_all="GammaZero"

# cau_cluster
module load cmake/3.18.4 
module load gcc/10.2.0 
module load eigen/3.3.8
module load gsl/2.5

# mkdir
rm -rf /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/0_bin
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/1_input
mkdir /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/2_output

# move files to rsync -v -a *.tar.gz sunbo481@nesh-fe.rz.uni-kiel.de:"/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/"

mv /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/2_RealData_mean.tar.gz /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/1_input 
cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/1_input/ && tar -xzf 2_RealData_mean.tar.gz


######################################
#       Dataset complete             #
######################################

#######################################################
#            all genes run for the mean               #
#######################################################

cd /gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/0_bin

# folders

input_FILE="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/1_input/2_RealData_mean"
output_FOLDER="/gxfs_work1/cau/sunbo481/projects/0_2020_POSTDOC_MAIN/$Species_is/3_analysis/16_grapes/9_recombination_15cat/3_singleRun/2_output"


ls $input_FILE/*.sfs | xargs -n 1 basename > SFS_run_for_mean.txt
sed -i 's/.sfs//' SFS_run_for_mean.txt
mkdir log

# create array batch file
touch jobs_SFS_run_for_mean.sh
while read line
do
echo "/gxfs_home/cau/sunbo481/software/BPP_danilo3/bin/grapes -in $input_FILE/$line.sfs -out $output_FOLDER/$line.res -model $Best_grapes_model_all -no_div_param" >> jobs_SFS_run_for_mean.sh
done < SFS_run_for_mean.txt

cat SFS_run_for_mean.txt | wc -l # this is the number of jobs in the array 300
max_array_num=$( cat SFS_run_for_mean.txt | wc -l )


# Create the submission script
echo -e '#!/bin/bash' > array_run_for_mean.sh
echo -e '# get sample' >> array_run_for_mean.sh
echo -e 'sample_ID=$(awk 'NR=='$SLURM_ARRAY_TASK_ID'' jobs_SFS_run_for_mean.sh)' >> array_run_for_mean.sh
echo -e '# run the command' >> array_run_for_mean.sh
echo -e '$sample_ID' >> array_run_for_mean.sh

echo "sbatch --job-name="$Abb_slurm" --array=1-$max_array_num --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=32G --error=log/log.%A_%a.err --output=log/log.%A_%a.out --partition=cluster array_run_for_mean.sh" > submit_1st.sh

bash submit_1st.sh
