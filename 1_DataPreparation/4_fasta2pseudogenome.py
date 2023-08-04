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
#               * standard gff annotation file (needs to be sorted)       #
#               * maffilter exon output file in fasta.                    #
#                                                                         #
###########################################################################
#
#
#
##################################### Script will do ##################################### 
# File input needed for pseudogenome preparation
# * STEP 1: Create a dataframe from GFF, with last column being gene_id repeated over exons (no output file).
# * STEP 2: Create a txt file with column 1 exon_id (contig_start-end) and column 2 gene_id as is (1 output file).
# * STEP 3: Create a txt tab-delimited file with column 1 exon_id (contig_start-end) and column 2 strand (+ or -) (1 output file).
# * STEP 4: Create a txt tab-delimited file with column 1 gene_id and column 2 strand (+ or -) (1 output file).
# * STEP 5: Create a file from maffilter exon list changing the header - Old header ">Ztritici-11(+)/1099-1510", new header ">11_1099-1510" (1 output file).
##########################################################################################

######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

####################################### STEP 1 ####################################### 
# This block will parse a GFF file into a GFF datafame containing only exons.
# Then will collect gene_id (name of gene e.g. Zt09_11_00001) for all exons (yes, repited gene_id's if more then 1 exon)
# The final goal is a dataframe GFF-like with gene_id repeated over exons (no output saved to file)
# IMPORTANT: The way the GFF is coded depends on software used, thus GFFs might be different and it can affect the script
# Example GFF file: /Users/danilo/Desktop/PostDoc/Stukenbrock/Project/2020_postdoc_main/03_results/20210216_mergeExon/01_run/8_RStudio_merge/Example_files/chr11.12.13.Ztritici_annotation_Grandaubert_noChr_noModel.sorted.gff
###########################################################################################
# load modules
import gffutils #http://daler.github.io/gffutils/index.html
import pandas as pd
from Bio import SeqIO
from Bio import Seq
import os
import glob
from collections import defaultdict
import csv


# Needs to be adjusted per species
Species_is="Verticillium"
GFF_name="Vdahliae_noDup_genomic.checked.gff"

# cd /home/pereira/2020_POSTDOC_MAIN/Verticillium/3_analysis/0_scripts/5_fasta2pseudogenome_CDS

## fixed
BASE=f"/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis"

# path to gff
GFF_file=f"{BASE}/4_maffilter/0_files_maffilter/{GFF_name}"

# Create a database. If already created, comment it out.
db = gffutils.create_db(GFF_file, dbfn=f'{GFF_file}.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

# Get GFF database
current_db = gffutils.FeatureDB(f'{GFF_file}.db', keep_order=True)

# Extract the exons from the database. This information is in the third collumn in the GFF.
CDS_list=[] # Creates an empty list, fill it with exons extracted from the gff file
for cds in db.features_of_type('CDS', order_by=('seqid','start')): # take 'exon' and order by chromosome then start position
    CDS_list.append(str(cds).split("\t")) # Collect all rows that are exon, and divide the row according to \t (space)

# take the EXON_list (rows collumns delimted by \t) and create dataframe with column names of a common gff
GFF_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
CDS_df = pd.DataFrame(CDS_list, columns=GFF_columns)

## IMPORTANT
## Up to this point a problem for scalating is related to how the GFF is coded. 
## Now the gene name is needed and it might also depend on how the GFF is coded.
# In this case of C. parasitica the last column for exons are like Parent=rna-XP_040780365.1, and the protein name is XP_040780365.1. What I can do is to erase `rna-`, which will leave Parent=XP_040780365.1 for exon's last column.
#CDS_df['attribute'] = CDS_df['attribute'].str.replace('rna-','',regex=True)
#CDS_df['attribute'] = CDS_df['attribute'].str.replace('[|,_]','',regex=True)

# get gene_id
CDS_df["gene_id"] = CDS_df["attribute"].str.split('=rna-', expand=True)[1]

#CDS_df
############## RC case, not sure if this is C parasitica case.
# I tried it, but then in STEP 3 of script 5, it will give an error because these tRNA exons will be missing. So I commented this part out, and proceeded without removing these at this point. Need to be removed later!
# OBS: there are some rna-xxx (tRNA exon). A total of 198 rows need to be removed. 3 containing pseudogenictranscript (e.g. pseudogenictranscript.SNOG401040A), and 195 containing rna-xxx.
#exon_df = exon_df[~exon_df.gene_id.str.contains("rna")]
#exon_df = exon_df[~exon_df.gene_id.str.contains("pseudogenictranscript")]
##############

#exon_geneID_df

####################################### STEP 2 ####################################### 
# Final goal is to output a txt tab-delimited file with column 1 exon_id (contig_start-end) and column 2 gene_id as is.
# Output saved to 1 file
########################################################################################### 

# Create a data frame with two columns. column 1 exon_id (contig_start-end) end and column 2 gene_id
CDS_geneID_df_temp = CDS_df[['seqname', 'start', 'end', 'gene_id']].copy()

## IMPORTANT
## Careful with 0-based or 1-based coordinates in GFF versus maffilter output. For more: https://www.biostars.org/p/84686/
# The position "end" given in the gff is one (1) base shorter then the end position output from maffilter (either exon fasta header or statistics file)
CDS_geneID_df_temp['end'] = CDS_geneID_df_temp['end'].astype(int) # ATTENTION to wrong conversions!! This will transform the values in column end (which are str) into "int", so addition is possible.
#exon_geneID_df_temp.dtypes
CDS_geneID_df_temp['end'] += 1 # add one base to the end positions
CDS_geneID_df_temp['end'] = CDS_geneID_df_temp['end'].astype(str) # convert the values back to string
#exon_geneID_df_temp

# create the exon_id vs gene_id key, to be used to change the headers in the fasta file to gene_id
CDS_geneID_df_temp["CDS_id"] = CDS_geneID_df_temp["seqname"] + "_" + CDS_geneID_df_temp["start"] + "-" + CDS_geneID_df_temp["end"] #It needs to match the maffilter output in the exon fasta sequence.
CDS_geneID_df = CDS_geneID_df_temp[['CDS_id', 'gene_id']]
#CDS_geneID_df


# export exon_id vs gene_id key file
CDS_geneID_df.to_csv(f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/CDS_geneID.txt",sep="\t", index=False, header=False)

####################################### STEP 3 ####################################### 
# Final goal is to output a txt tab-delimited file with column 1 exon_id (contig_start-end) and column 2 strand (+ or -).
# Output saved to 1 file
########################################################################################### 

strand_CDSid_df = CDS_df[['seqname', 'start', 'end', 'strand']].copy() # select columns to be used from exon_df

## IMPORTANT
## Careful with 0-based or 1-based coordinates in GFF versus maffilter output. For more: https://www.biostars.org/p/84686/
# The position "end" given in the gff is one (1) base shorter then the end position output from maffilter (either exon fasta header or statistics file)
strand_CDSid_df['end'] = strand_CDSid_df['end'].astype(int) # ATTENTION to wrong conversions!! This will transform the values in column end (which are str) into "int", so addition is possible.
#exon_geneID_df_temp.dtypes
strand_CDSid_df['end'] += 1 # add one base to the end positions
strand_CDSid_df['end'] = strand_CDSid_df['end'].astype(str) # convert the values back to string
#strand_exonid_df

# build exonid again.
strand_CDSid_df["CDS_id"] = strand_CDSid_df["seqname"] + "_" + strand_CDSid_df["start"] + "-" + strand_CDSid_df["end"]
strand_CDSid_df = strand_CDSid_df[['CDS_id', 'strand']]
#strand_exonid_df

# export new file
strand_CDSid_df.to_csv(f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/CDS_strand.txt",sep="\t", index=False, header=False)


####################################### STEP 4 ####################################### 
# Final goal is to output a txt tab-delimited file with column 1 geneID and column 2 strand (+ or -).
# Output saved to 1 file
########################################################################################### 

# copy desired columns to new DF
strand_geneID_df = CDS_df[['gene_id', 'strand']].copy()
strand_geneID_df

# collapse geneID to unique names
strand_geneID_df = strand_geneID_df.drop_duplicates(subset=['gene_id'])
strand_geneID_df

# export new file
strand_geneID_df.to_csv(f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/geneID_strand.txt",sep="\t", index=False, header=False)

# set directory with maffilter CDS sequences. OBS: one CDS fasta per isolate
Old_Header_Directory = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/1_maffilter_output/'
new_Header_Directory = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/2_fasta2pseudo_newheader/'


# start loop to catch the path of input fasta files
for fasta_file in glob.glob(os.path.join(Old_Header_Directory, '*.fasta')):
    #print(fasta_file)
    #print(os.path.basename(fasta_file))
    sample_fasta = os.path.basename(fasta_file)
    sample_ID = sample_fasta.split(".")[0]
    output_file = new_Header_Directory + sample_ID + ".filtered.newheader.fasta"
    #name_modify = "desktop/file/danilo/the_file_is:%s" % sample_fasta
    #output_path = new_Header_Directory + "/" + 
    print(fasta_file)
    sequencias_new_header = [] # create an empty list
    # this loop will open each fasta, and change the sequences header
    for sequencia in SeqIO.parse(fasta_file, 'fasta'): # parse the fasta using SeqIO
        chromosome_seqname = str(sequencia.id).split("/")[0][:-3].split("-")[1] # Get the sequence id, split, take first column (0), remove last 3 characters ([:-3]), split again based on -
        #print(chromosome_seqname)
        start_end_list = str(sequencia.id).split("/")[1] # will take '1099-1510'
        #print(start_end_list)
        new_header_fasta = chromosome_seqname + "_" + start_end_list # merge to compose the header
        #print(header_fasta)
        sequencia.id = new_header_fasta # assign the new header to the id
        sequencia.description = "" # this strips the old header out, otherwise old + new will be merged
        sequencias_new_header.append(sequencia) # append the current sequence to a list
    SeqIO.write(sequencias_new_header, output_file, 'fasta')

##############################################################################################################################

####################################### Script will do #################################### 
# Generate pseudogenome
# * STEP 1: Create a series of dictionaries needed to reverse complement exons on negative strand
# * STEP 2: Will perform the reverse complementation of exons on negative strand
# * STEP 3: Will merge exons with equal name, rename the exons to gene name, and output one fasta per isolate with genes in it
###########################################################################################

#################################### Necessary files ######################################
# * A txt tab-delimited file with column 1 CDS_id (contig_start-end) and column 2 gene_id as is (1 output file) (from Step 2, script1)
# * A txt tab-delimited file with column 1 CDS_id (contig_start-end) and column 2 strand (+ or -) (1 output file) (from Step 3, script1)
# * A txt tab-delimited file with column 1 gene_ID and column 2 strand (+ or -) (1 output file) (from Step 4, script1)
# * A fasta file containing exons, headers as CDS_id. (from Step 5, script1) ** moved to script 1.1
###########################################################################################

# load all necessary key files
CDS_geneID_key = f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/CDS_geneID.txt"
CDS_strandID_key = f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/CDS_strand.txt"
GENEID_strandID_key = f"{BASE}/4_maffilter/0_files_maffilter/1_key_files/geneID_strand.txt"
print("Key files loaded")

# folder containing maffilter exon output with new header (header as exonID)
#new_Header_Directory = "/home/pereira/2020_POSTDOC_MAIN/Penicillium/3_analysis/4_maffilter/2_fasta2pseudo_newheader"

# output folder for non-final alignment
output_temp = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/3_temp/'

# final output folder
output_final_seq = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/4_fasta2pseudo_CDS/'

print("Path to folders loaded")

####################################### STEP 1 ####################################### 
# This block will take the "exon:gene" key and create a dictionary (no output file).
# This block will take the "exon:strand" key and create a dictionary (no output file).
# This block will take the "geneID:strand" key and create a dictionary (no output file).
# Final goal is to produce dictionaries needed to map exons to genes and strand.
####################################### 

# "exonID:geneID" dictionary
# This block will open the CDS_geneID_key file, and create a tuple (list not modifiable) dictionary "exonID:geneID"
with open(CDS_geneID_key) as f: # standard opening of a file in python. Read.
    #print(f)
    lines = f.read().splitlines()
    #print(lines)
    CDS_to_gene = dict(tuple(line.split()) for line in lines) # since there are two values per line, create a dict for every line in list (lines)

# "exonID:strandID" dictionary
with open(CDS_strandID_key) as f:
    #print(f)
    lines = f.read().splitlines()
    #print(lines)
    CDS_to_strand = dict(tuple(line.split()) for line in lines)

# "geneID:strandID" dictionary
with open(GENEID_strandID_key) as f:
    #print(f)
    lines = f.read().splitlines()
    #print(lines)
    geneID_to_strand = dict(tuple(line.split()) for line in lines)

print("Tuple dictionaries created")

####################################### STEP 2 ####################################### 
# This block will take the exons output from maffilter and reverse complemente those on negative strand
# Final goal is to reverse complement exons on negative strand.
###################################################################################### 
print("Starting reverse complementation")
# start loop
for fasta_file in glob.glob(os.path.join(new_Header_Directory, '*.fasta')): # take all fasta files in folder
    #print(fasta_file)
    #print(os.path.basename(fasta_file))
    sample_fasta = os.path.basename(fasta_file) # take only the part after last /
    #print(sample_fasta)
    print("Reverse complementing CDS from sample: ", sample_fasta)
    CDS_maffilter = fasta_file


    # will open CDS_maffilter fasta, extract id, compare id with dictionary, and if negative reverse_complement and replace sequence, then append sequence. If not == "-" just append sequence.
    new_CDS_record = [] # create empty list the append sequence
    for sequencia in SeqIO.parse(CDS_maffilter, 'fasta'):
        for key, value in CDS_to_strand.items():
            if key == sequencia.id and value == "-":
                sequencia = sequencia.reverse_complement(id=True, name=True, description=True, annotations=False)
        new_CDS_record.append(sequencia)

    #print(new_CDS_record)

    # export exon sequence as fasta for future checking if necessary and openning in Step 3
    # save the new list (new_CDS_record) as fasta file
    sample_ID = sample_fasta.split(".")[0] # extract only isolate ID (first part before first .)
    output_file = output_temp + sample_ID + ".filtered.newheader.rc.fasta" # construct whole path and file name for output.

    with open(output_file, 'w') as g:
        SeqIO.write(new_CDS_record, g, 'fasta')

print("Reverse complementation finished")

####################################### STEP 3 ####################################### 
# This block will replace exonID by geneID, having geneId as many times there are exons
# Then will merge everything with a same geneID
# Final goal is to output a fasta per isolate with complete merged list of exons
############################################################################ 
print("Starting merge of CDS")

# start loop
for fasta_file in glob.glob(os.path.join(output_temp, '*.fasta')): # take all fasta files in folder
    print(fasta_file)
    #print(os.path.basename(fasta_file))
    sample_fasta = os.path.basename(fasta_file) # take only the part after last /
    #print(sample_fasta)
    print("Replacing cdsID by geneID and merging CDSs for sample: ", sample_fasta)
    new_CDS_record = fasta_file

    # this line will open the previously reverse complemented fasta and create a dictionry (header, sequence)
    new_CDS_records = SeqIO.to_dict(SeqIO.parse(new_CDS_record, 'fasta'))
    #print(new_CDS_records)
    
    # replace CDS_id (current header) by gene_id
    # this loop will take the CDS_to_gene file (exonid, geneid) and use that to replace the header of new_CDS_records fasta file, which now contains CDS_id, by the geneid.
    gene_records = defaultdict(list) # explanation on defaultdict() https://stackoverflow.com/questions/5900578/how-does-collections-defaultdict-work
    for CDS_id, record in new_CDS_records.items():
        gene_id = CDS_to_gene[CDS_id]
        #print(gene_id)
        gene_records[gene_id].append(record)

    #print(gene_records)
    
    # merge exons (everything with same gene_id)
    # this loop will merge equal geneids, which represets same exons of a gene
    cds_records = []
    for gene_id in gene_records:
        gene_record = gene_records[gene_id][0]
        if len(gene_records[gene_id]) > 1:
            sequence = ''.join([str(gene_record.seq) for gene_record in gene_records[gene_id]])
            gene_record.seq = Seq.Seq(sequence)
        gene_record.id = gene_id
        gene_record.description = gene_id # see https://www.biostars.org/p/156428/
        cds_records.append(gene_record)
    #print(cds_records)

    # will open previous cds_records, extract id, compare id with dictionary, and if negative reverse_complement and replace sequence, then append sequence. If not == "-" just append sequence.
    new_cds_records = [] # create empty list the append sequence
    for sequencia in cds_records:
        for key, value in geneID_to_strand.items():
            if key == sequencia.id and value == "-":
                sequencia = sequencia.reverse_complement(id=True, name=True, description=True, annotations=False)
        new_cds_records.append(sequencia)

    # export merged exons
    sample_ID = sample_fasta.split(".")[0] # extract only isolate ID (first part before first .)
    output_file = output_final_seq + sample_ID + ".completeCDS.fasta" # construct whole path and file name for output.

    with open(output_file, 'w') as g:
        SeqIO.write(new_cds_records, g, 'fasta')

print("Merge of CDS finished")

# split fasta into sequences, with new seq. header as sample, and append matching geneID. https://www.biostars.org/p/273248/
# OBS: output will be written to the current directory (not necessarily the one where script is)


####################################################
#### This script needs to run for all samples   ####
####################################################

####################################### Script will do #################################### 
# * STEP 1: Take individual fasta per isolate, output as many fasta files as genes, with individuals inside
###########################################################################################

#################################### Necessary files ######################################
# * Individual fasta per isolate, with geneID header (merged exon already)
###########################################################################################

####################################### STEP 1 #######################################
# Will take all fasta in a folder
# Will extract from each fasta the header (id) and sequence
# Will writte/append sequences to fasta named after the geneID, with new sequence headers as sampleID
# If one isolate is missing a gene, the gene fasta will not have that gene and loop will not crash
########################################################################################### 

# folder containing the fasta files per sample
input_Directory = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/4_fasta2pseudo_CDS/'

# folder for output
output_folder = f'/home/pereira/2020_POSTDOC_MAIN/{Species_is}/3_analysis/4_maffilter/6_split_into_genes_fasta.v2'

print("Folders loaded")

print("Loop starting")
# loop start
for fasta_file in glob.glob(os.path.join(input_Directory, '*.fasta')):
    #print(fasta_file)
    sample_fasta = os.path.basename(fasta_file)
    sample_ID = sample_fasta.split(".")[0]
    print("Splitting sample: ", sample_fasta)

    f_open = open(fasta_file, newline='\n')

    for rec in SeqIO.parse(f_open, "fasta"):
        output_file = os.path.join(output_folder, rec.id + ".fasta") # create the outgroup naming
        seq = rec.seq
        id_file = open(output_file, "a+") # a+ > append
        id_file.write(">"+str(sample_ID)+"\n"+str(seq)+"\n") # \n > new line
        id_file.close()


    f_open.close()

print("Splitting sequences finished. Bye Danilo!")
print("Script done. Bye Danilo!")



######################################################################################################################
#                                                      END                                                           #
######################################################################################################################