
__author__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"
# coding: utf-8

import sys
import csv
from Bio import Entrez
from collections import Counter
import matplotlib.pyplot as plt

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
    with open(output_file, mode='w') as efetch_output_file:
        efetch_output = csv.writer(efetch_output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        efetch_output.writerow(['ID', 'Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name', 'Organism', ' Strain', 'Assembly'])

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

            organism = fields[8]  # assuming Organism is the 9th field

            # Extract species name from organism field
            species = organism.split(" ")[0]

            organisms[species] += 1  # increment the count for this species

        return organisms

def main(input_file, output_file, organisms_file="organisms.txt"):
    input_handle = fetch_data(input_file)
    organisms = process_data(input_handle, output_file)

    print("Organisms by count:")
    for organism, count in organisms.most_common():
        print(f"{organism}: {count}")

    with open(organisms_file, "w") as f:
        for organism, count in organisms.most_common():
            f.write(f"{organism}: {count}\n")

    # Filter organisms with count > 1
    filtered_organisms = Counter({organism: count for organism, count in organisms.items() if count > 2})

    # Generate and save bar chart with higher resolution and adjusted size
    plt.figure(figsize=(16, 12), dpi=300)  # Adjust the figsize and dpi as needed
    organisms_names = [organism for organism, _ in filtered_organisms.most_common()]  # Use filtered organisms
    organisms_counts = [count for _, count in filtered_organisms.most_common()]  # Use filtered organisms
    plt.barh(organisms_names, organisms_counts)
    plt.xlabel("Count")
    plt.ylabel("Genus")
    plt.title("Distribution of Bacteria with Count > 1")
    plt.tight_layout()  # Adjust spacing between plot elements
    plt.savefig("bacteria_distribution.png")

    print("Program finished")

if __name__ == "__main__":
    main(sys.argv[1], "efetch_output.tsv")
