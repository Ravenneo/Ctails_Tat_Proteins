#!/usr/bin/env python
# coding: utf-8


__author__ = 'Dr Giusy Mariano'
__email__ = 'giusy.mariano@ncl.ac.uk'
__license__ = "GPL"
__modification__ = 'Jose Jesus Gallego-Parrilla'
__email__ = 'J.J.Gallego-Parrilla2@newcastle.ac.uk'
__license__ = "GPL"

# Import necessary modules.
# sys module is used to work with command-line arguments.
# pandas is used for handling data in structured format (dataframes).
import sys
import pandas as pd

# The first command-line argument is the name of the input file.
file_name = sys.argv[1]

# The second command-line argument is the base name of the output files.
file_name_output = sys.argv[2]

# Read the input file into a pandas dataframe.
# low_memory=False is used to eliminate a warning that can occur when inferring data types.
# names provides column names for the dataframe.
df = pd.read_csv(file_name, sep="\t", low_memory=False, 
                 names=['ID', 'Source', 'Nucleotide Accession', 'Protein', 'Protein Name', 'Start', 
                        'Stop', 'Strand', 'Organism',' Strain', 'Assembly'])

# Get names of indexes for which rows have to be dropped
# Rows where the 'Protein Name' is 'hypothetical protein' will be dropped.
indexNames = df[df['Protein Name'] == 'hypothetical protein'].index

# Delete these row indexes from dataFrame
df.drop(indexNames, inplace=True)

# Write out the dataframe to a tsv file.
df.to_csv(file_name_output + ".tsv", sep="\t", index=False)

# Create a new dataframe where only one representative for each 'Protein' is kept.
# Rows are sorted by 'Protein', and for each unique 'Protein', only the first occurrence is kept.
df2 = df.sort_values(by="Protein", axis=0, ascending=True, inplace=False).drop_duplicates(subset=['Protein'],keep='first')

# Write out this new dataframe to a tsv file.
df2.to_csv(file_name_output + "_one_WP_per_assembly.tsv", sep="\t", index=False)

print ('program finished')






