################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
###########################################################################
# This script will prepare all necessary files for pseudogenome assembly  #
#                                                                         #
# Files needed:                                                           #
#               * protein fasta file from target and outgroup species     #
#                                                                         #
###########################################################################
#
#
#
############################################## Script will do ################################################## 
# File input needed for pseudogenome preparation                                                               #
# * this script read a single fasta file with multiple sequences and output a single file per protein sequence #
#                                                                                                              #
################################################################################################################              

######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

# load module
import pandas as pd
from Bio import SeqIO
import glob
import os

# set names per species
Species_is="Verticillium"
Abb_is="Vdahliae"

print("Importing data")
# import results from orthofinder
SCO_key = f"/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/8_select_orthologs/OrthoFinder/Results_Mar10/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
matching_genes_raw = f"/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/8_select_orthologs/OrthoFinder/Results_Mar10/Orthogroups/Orthogroups.txt"

# load data
SCO_key_list = (pd.read_table(SCO_key, header = None))
#SCO_key_list.rename(columns={ SCO_key_list.columns[0]: "orthogroup" }, inplace = True)

matching_genes_raw_DF = (pd.read_table(matching_genes_raw, header = None, sep=" ", engine='python'))
#matching_genes_raw_DF.rename(columns={ matching_genes_raw_DF.columns[0]: "orthogroup" }, inplace = True)
matching_genes_raw_DF[0] = matching_genes_raw_DF[0].str.replace(r':', '')

# merge
# excellent: https://stackoverflow.com/questions/53645882/pandas-merging-101
SCO_DF = SCO_key_list.merge(matching_genes_raw_DF, on=0, how='left')
SCO_DF = SCO_DF[[0,1,2]]

SCO_DF_key = SCO_DF[[1,2]]

print("Exporting")

# export data
SCO_DF.to_csv(f"/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/8_select_orthologs/Species_OUT_orthologs.txt", sep='\t', index=False, header=None)
SCO_DF_key.to_csv(f"/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/8_select_orthologs/Species_OUT_orthologs_keyOnly.txt", sep='\t', index=False, header=None)

print("Splitting outgroup")

# folder for output
output_folder = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/9_outgroup_genes'

# multisequence fasta file
single_fasta = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/8_select_orthologs/augustus.hints.codingseq.fasta'

# this loop will read trough the multisquence fasta file and extract sequences into files
f_open = open(single_fasta, newline='\n')
for rec in SeqIO.parse(f_open, "fasta"):
    print("Splitting CDS: ", rec.id)
    output_file = os.path.join(output_folder, rec.id + ".fasta") # create the outgroup naming
    seq = rec.seq
    id_file = open(output_file, "a+") # a+ > append
    id_file.write(">"+str(rec.id)+"\n"+str(seq)+"\n") # \n > new line
    id_file.close()

print("Splitting sequences finished. Bye Danilo!")

print("Bye Danilo!")


######################################################################################################################
#                                                      END                                                           #
######################################################################################################################