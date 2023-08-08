################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
#############################################################################################################
# This script is intended to itinerate over the gap scores and fail alignments given a threshold            #
#                                                                                                           #
# Files needed:                                                                                             #
#             * Output from script 10_Gap_threshold.v2.sh                                                   #
#                                                                                                           #
#                                                                                                           #
#############################################################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################


# this python script will parse the dataframe generated in the previous script and give file names that fail a gap treshold

# load libs
import pandas as pd
import numpy as np

# give here the threshold as e.g. 5 (means 5% gap)
gap_threshold=5

# load data. FILE NEEDS TO BE IN SAME FOLDER AS EXECUTION
gap_DF = pd.read_csv("GAP_quantification.txt", sep='|', header=None)
#gap_DF.columns

# subset and change column names
gap_simple_DF = gap_DF[[0, 1, 3, 5,6,7,8, 9, 13]] # select columns to subset
gap_simple_DF.columns=["file","seq", "type", "num_seq", "sum_len", "min_len", "avg_len", "max_len", "sum_gap"] # rename columns
gap_simple_DF.columns

# select rows that contain the string FASTA in the column type
gap_simple_DF = gap_simple_DF[gap_simple_DF['type'].str.match('FASTA')] 
gap_simple_DF

# change datatype, transfor data to numeric
gap_simple_DF[["num_seq", "sum_len", "sum_gap"]] = gap_simple_DF[["num_seq", "sum_len", "sum_gap"]].apply(pd.to_numeric)
gap_simple_DF.dtypes

# this code uses numpy to create a new column called gap_percentage
# this new column will be the number in sum_len if sum_gap is < 1, or if > 1 will be ( sum_gap / sum_len ) * 100
# meaning, if no gap in sequence, copy sum_gap value, otherwise, do ( sum_gap / sum_len ) * 100 
gap_simple_DF['gap_percentage'] = np.where(gap_simple_DF['sum_gap'] < 1, gap_simple_DF['sum_gap'], (gap_simple_DF['sum_gap']/gap_simple_DF['sum_len'])*100)
gap_simple_DF

# create a list of our conditions. Basically if smaller or equan then the threshold, print PASS to the column, otherwise, print FAIL
conditions = [
    (gap_simple_DF['gap_percentage'] <= gap_threshold),
    (gap_simple_DF['gap_percentage'] > gap_threshold)
    ]

# create a list of the values we want to assign for each condition
values = ['PASS', 'FAIL']

# create a new column and use np.select to assign values to it using our lists as arguments
gap_simple_DF['FLAG'] = np.select(conditions, values)

# select file names that failed
fail_gap_DF = gap_simple_DF[gap_simple_DF['FLAG'].str.match('FAIL')] 
fail_gap_DF = fail_gap_DF[["file", "FLAG"]] # select columns to subset

# select file names that failed and get unique names only
gap_simple_uniqueONLY_DF = fail_gap_DF.drop_duplicates('file')
gap_simple_uniqueONLY_DF

# output analysis

# full output
gap_simple_DF.to_csv("1_GAP_summary.txt", sep='\t',index=False)

# file names which failed only
fail_gap_DF.to_csv("2_GAP_geneNAME_failed.5.txt", sep='\t',index=False,header=False)

# file names which failed and unique only
gap_simple_uniqueONLY_DF.to_csv("3_GAP_geneNAME_failed.5.UNIQUEnames.txt", sep='\t',index=False,header=False)


print("Bye Danilo")

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
