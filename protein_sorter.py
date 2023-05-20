#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses the ProteinHtml class to find all unique protein names in an HTML file and
write them into a CSV file along with their counts of occurrences.

__author__ = 'Jorge Camarero Vera, Jose Jesus Gallego-Parrilla'
__license__ = "GPL"
__maintainer__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"
"""

# Import the pandas library for data manipulation
import pandas as pd

# Import the ProteinHtml class from the ProteinHtml module
from ProteinHtml import ProteinHtml 

# Create a new instance of the ProteinHtml class with the "chandra.html" file
protein = ProteinHtml("chandra.html")

# Find all unique protein names in the file
found = protein.findAllProteinNames()

# Create two lists to hold the protein names and their counts of occurrences
names = []
places = []

# Iterate over all proteins found in the file
for elem in found:
  names.append(elem)  # Append the protein name to the names list
  places.append(found[elem])  # Append the count of occurrences to the places list

# Create a dictionary to hold the names and places lists
data = {'Name' : names, 'Place' : places}

# Create a DataFrame from the dictionary
dataframe = pd.DataFrame(data)
df = pd.DataFrame(data)

# Write the DataFrame to a CSV file
df.to_csv('all_proteins.csv', index=False)

# Print the DataFrame
print(dataframe)
