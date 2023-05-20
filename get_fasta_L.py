#!/usr/bin/env python
# coding: utf-8


__author__ = 'Jose Jesus Gallego-Parrilla'
__license__ = "GPL"
__maintainer__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"


# Import necessary modules.
# sys module is used to work with command-line arguments.
# csv module is used to read from and write to csv files.
# Entrez from Bio is used to fetch data from NCBI databases.
import sys
import csv
from Bio import Entrez

# Set the email address to be used by NCBI when you're accessing their database.
Entrez.email = "your@email.com"

# Check if the correct number of command-line arguments has been provided.
# If not, print a usage message and exit.
if len(sys.argv) < 3:
    print("It needs input and output arguments.")
    print("Ej: python get_fasta_L.py efetch.csv efetch_output.txt")
    sys.exit(0)

# The first command-line argument is the name of the csv file containing the accession numbers.
csv_file = sys.argv[1]

# The second command-line argument is the name of the output text file where the fetched sequences will be written.
txt_file = sys.argv[2]

# Initialize an empty list to store the accession numbers.
list_of_accession = []

# Open the csv file and read the accession numbers, appending them to the list.
with open(csv_file, 'r', encoding='utf-8-sig') as csvfile:
    efetchin = csv.reader(csvfile, delimiter=',')
    for row in efetchin:
        list_of_accession.append(str(row[0]))

# Open the output text file in write mode.
with open(txt_file, mode='w') as efetch_output:

    # Fetch the protein sequences in FASTA format from NCBI database.
    input_handle = Entrez.efetch(db="protein", id=list_of_accession, rettype="fasta")

    # Open the output text file in append mode.
    output_handle = open(txt_file, "a")

    # Write the fetched sequences to the output file.
    for line in input_handle:
        output_handle.write(line)

# Close the handles.
input_handle.close()
output_handle.close()

print('program finished')
