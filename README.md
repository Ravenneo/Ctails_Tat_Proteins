![Screenshot from 2023-05-20 21-22-16](https://github.com/Ravenneo/Ctails_Tat_Proteins/assets/41577767/eee0e72a-be33-4729-b6df-f3c1fe61a2c2)


# Ctails_Tat_Proteins
A repository containing the scripts used in my thesis: 

**Integration of tail-anchored proteins by the twin-arginine translocas**


This repository contains Python scripts for bioinformatics analysis used in my thesis, specifically designed for protein sequence and assembly data. The scripts were authored and modified by Jose Jesus Gallego-Parrillam Dr. Giusy Mariano and Jorge Camarero Vera. The scripts were used wiht the data provided by Profesor Tracy Palmer from Palmer Lab (Newcastle University)

Content

    ·cala.py: Extracts protein IDs from "chandra.html" for console viewing or CSV export.
    
    ·get:_fasta_L.py: Fetches protein sequences from the NCBI database given a CSV list of accession numbers, and writes them to a text file.
    
    ·reptile2.py: Processes a tab-separated file to remove 'hypothetical protein' entries and keeps only the first occurrence of each unique protein.
    ·Main_proteins.py: Uses the ProteinHtml class to extract protein-related information from an HTML file.
    
    ·ProteinHtml.py: Defines ProteinHtml class for extracting protein data from HTML files.
    
    ·Protein_sorter.py: Creates a CSV file of unique protein names and their counts using the ProteinHtml class.
    
    ·Incomplete_protein_counter.py: Finds and counts occurrences of specific protein names/terms in an HTML file.
    
    ·Arcana.py: Searches for Hidden Markov Models on protein sequences and generates a CSV of hits. Also fetches protein details from NCBI's Entrez service.
    
    ·Etna2.py: Standalone version of Arcana's second part. Fetches protein information from NCBI's Entrez service given a list of accessions, and outputs a detailed table.


These scripts are licensed under the GNU General Public License (GPL). For more information, please see the LICENSE file in this repository.


# UPDATE
20-July-2023

## wanda.py 
This Python script is designed for the purpose of retrieving and analyzing protein data from the NCBI databases. The script uses the Biopython library's Entrez module to access NCBI's Entrez databases, and fetches data about proteins based on their accession numbers.

The script takes as input a CSV file containing a list of protein accession numbers. It uses the Entrez efetch method to fetch the data in TSV format for each protein in the list.

The fetched data includes the following fields: 'ID', 'Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name', 'Organism', 'Strain', and 'Assembly'. The script writes this data to a new TSV file named efetch_output.tsv.

During the data fetching and writing process, the script applies two filters to the data:

    It ignores any proteins that are labeled as "hypothetical proteins".
    It ensures that each protein ID only appears once in the output file, ignoring subsequent occurrences of the same ID.

The script also keeps track of the organisms associated with the proteins. It counts the occurrences of each organism and writes these counts to a new file named organisms.txt. In both the console output and the file, the organisms are listed in descending order of their counts, from most to least.

Finally, the script prints a message to the console to indicate that it has finished running.

The script requires the Biopython and pandas libraries to run. It should be run from the command line with the name of the input CSV file as the argument. The email address used for the Entrez database queries is hard-coded into the script as "raven.neo@gmail.com" and may need to be changed to the user's own email address.

Example usage: python3 script.py input.csv
