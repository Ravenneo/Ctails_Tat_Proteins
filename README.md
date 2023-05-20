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
