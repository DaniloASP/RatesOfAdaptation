################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
########################################################################
#   This script will take the output of LDhat and parse it             #
#                                                                      #
# Files needed:                                                        #
#             * Output files from LDhat parsed in the previous script  #
#                                                                      #
########################################################################
#
#
#
#####################################################################

# load libs
import pandas as pd
import numpy as np
import glob
import os

species_name="Zymoseptoria"
output_folder = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/13_Ldhat/7_summ"

#contig_info_folder = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/0_scripts/15_LDhat/0_first_script"
#input_files = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/13_Ldhat/6_stat"

#contig_length_file=f"{contig_info_folder}/contig_list_length.txt"

window_size_bp = int(20000)

# after summarizing all raw runs of LDhat, create the average from all runs and output a single file for the species
appended_data = []
for file in (glob.glob(f"{output_folder}/summ.*")):
  #print(file)
  rho_summ_df = pd.read_csv(file, sep="\t", header=None)
  rho_summ_df = rho_summ_df.tail(rho_summ_df.shape[0] -1) # remove header row
  # store DataFrame in list
  appended_data.append(rho_summ_df)

appended_data = pd.concat(appended_data)
appended_data = appended_data.rename(columns={0: 'window', 1: 'rho_mean', 2: 'contig'})

# summarize the data
appended_data["rho_mean"] = pd.to_numeric(appended_data["rho_mean"])
#appended_data["rho_mean"] = appended_data["rho_mean"].round(5) # round to 5 decimals
#appended_data["rho_mean"] = np.round(appended_data["rho_mean"], decimals=5)
#appended_data.dtypes

final_summ_data_DF = pd.DataFrame(appended_data.groupby(['window','contig'])['rho_mean'].mean())
final_summ_data_DF.reset_index(inplace=True)
#final_summ_data_DF
# export file
final_summ_data_DF.to_csv(f"{output_folder}/Final.{species_name}.LDhat.20kb.txt", sep='\t', index=False)


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

