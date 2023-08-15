################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################
#   This script will take the output of LDhat and parse it         #
#                                                                  #
# Files needed:                                                    #
#             * Output files from LDhat                            #
#                                                                  #
####################################################################
#
#
#
#####################################################################

# load libs
import pandas as pd
import numpy as np
import glob
import os

contig_info_folder = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/0_scripts/15_LDhat/0_first_script"
input_files = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/13_Ldhat/6_stat"
output_folder = "/home/pereira/2020_POSTDOC_MAIN/Zymoseptoria_RUN2/3_analysis/13_Ldhat/7_summ"

contig_length_file=f"{contig_info_folder}/contig_list_length.txt"

window_size_bp = int(20000)

# loop over LDhat output files and proceed with analysis for each
for file in (glob.glob(f"{input_files}/rep_*/*")):
  file_name = os.path.basename(file)
  species_name = file_name.split(".")[1]
  contig_code = int(file_name.split(".")[2]) # MAKE SURE THE CONTIG CODE IS A INT AND NOT A STRING!!!
  rep_code = file_name.split(".")[0]
  print(file_name)
  # create table with windows of 20kb
  contig_length_df = pd.read_csv(contig_length_file, sep="\t", header=None)
  contig_length_dict = contig_length_df.set_index(0)[1].to_dict() # convert to dictionary
  # take information from each contig
  for key, value in contig_length_dict.items(): # key is the contig and value the contig length
    if key == contig_code: # if current contig matches the contig in the dictionary, take the length
      #print(value)
      ### create a range df
      start = list(range(0,value,window_size_bp))
      end = start.copy()
      end[:] = [number - 1 for number in end]
      end.pop(0) # remove first element
      start = start[:-1] # remove last element
      range_chr= {'start':start, 'end':end}
      range_chr_df = pd.DataFrame(range_chr)
      print(file_name)
      # load file with LDhat output
      rho_raw = pd.read_csv(file, sep="\t")
      rho_raw = rho_raw.tail(rho_raw.shape[0] -1)
      rho_raw['window_range'] = 'range' # create a new column filled with 'range'
      for row_ID, row_content in rho_raw.iterrows():
        for row_ID2, row_content2 in range_chr_df.iterrows():
          if row_content[0] >= row_content2[0] and row_content[0] <= row_content2[1]: # row_content[0] is the locus, row_content2[0] the     start position and row_content2[1] the end position
          #range_to_append.append((row_content2[0],row_content2[1]))
            current_range = '-'.join(str(e) for e in (row_content2[0],row_content2[1])) # will convert both coordinates (that are a list    ), into a string, separated by '-', so it can be add to the row given by row_ID in the nex line
            rho_raw.loc[row_ID,'window_range'] = current_range # add current_range to the df
          #print("pass")
    
      # summarize the data
      summ_data_DF = pd.DataFrame(rho_raw.groupby("window_range")['Mean_rho'].mean())
      summ_data_DF.reset_index(inplace=True)
      summ_data_DF
  
      # create another loop to add coordinates where rho calculation was not possible
      for row_ID2, row_content2 in range_chr_df.iterrows():
        current_range = '-'.join(str(e) for e in (row_content2[0],row_content2[1]))
        #print(current_range)
        if summ_data_DF['window_range'].str.fullmatch(current_range).any(): # check if range is in the df, if yes will pass, if not, create a new row to be added
          pass
        else:
          new_line = {'Mean_rho': 0, 'window_range': (current_range)} # compose a new two column row
          summ_data_DF = summ_data_DF.append(new_line, ignore_index=True)   
      
      # add chromosome column
      summ_data_DF["contig"] = contig_code
      
      # export file
      summ_data_DF.to_csv(f"{output_folder}/summ.{file_name}", sep='\t', index=False)

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################




