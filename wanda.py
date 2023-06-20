
__author__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"
# coding: utf-8

import sys
import csv
from Bio import Entrez
from collections import Counter

def fetch_data(input_file, email="raven.neo@gmail.com"):
    Entrez.email = email
    list_of_accession = []
    with open (input_file, 'r') as csvfile:
        efetchin=csv.reader(csvfile, delimiter = ',')
        for row in efetchin:
            list_of_accession.append(str(row[0]))
        
    input_handle = Entrez.efetch(db="protein", id= list_of_accession, rettype="ipg", retmode="tsv")
    return input_handle

def process_data(input_handle, output_file="efetch_output.tsv"):
    with open(output_file, mode = 'w') as efetch_output_file:
        efetch_output = csv.writer(efetch_output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        efetch_output.writerow(['ID','Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name', 'Organism', ' Strain', 'Assembly'])
        
        organisms = Counter()  # use Counter instead of set
        seen_ids = set()  # to keep track of seen IDs
        for line in input_handle:
            decoded_line = line.decode("utf-8")  # decode bytes to string
            
            # split the line into fields
            fields = decoded_line.split("\t")
            
            # check if this protein is hypothetical
            protein_name = fields[7]  # assuming Protein Name is the 8th field
            if "hypothetical protein" in protein_name:
                continue  # skip this line
            
            # check if this ID has already been seen
            id_ = fields[0]  # assuming ID is the 1st field
            if id_ in seen_ids:
                continue  # skip this line
            
            seen_ids.add(id_)
            efetch_output_file.write(decoded_line)
            organisms[fields[8]] += 1  # increment the count for this organism

        return organisms

def main(input_file, output_file, organisms_file="organisms.txt"):
    input_handle = fetch_data(input_file)
    organisms = process_data(input_handle, output_file)

    print("Organisms by count:")
    for organism, count in organisms.most_common():  # use most_common() to order by count
        print(f"{organism}: {count}")  # print organism and count

    with open(organisms_file, "w") as f:
        for organism, count in organisms.most_common():
            f.write(f"{organism}: {count}\n")

    print ('program finished')

if __name__ == "__main__":
    main(sys.argv[1], "efetch_output.tsv")
