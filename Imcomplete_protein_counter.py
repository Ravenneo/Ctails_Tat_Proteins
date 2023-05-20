#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads an HTML file and counts the occurrences of specific proteins or protein-related terms.

__author__ = 'Jorge Camarero Vera, Jose Jesus Gallego-Parrilla'
__license__ = "GPL"
__maintainer__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"
"""

# Define a function to count the occurrences of a specific word in a file
def countWords(address, word):
    # Open the file in read mode
    logfile = open(address, "r")
    
    # Initialize the word count to 0
    wordcount = 0
    
    # Iterate over each line in the file
    for lines in logfile:
        # If the word is in the line, increment the word count
        if word in lines.split():
            wordcount += 1
    
    # Return the final word count
    return wordcount

# Import the regular expressions (re) module
import re

# Open the "chandra.html" file and read its contents
with open("chandra.html") as f:
    contents = f.read()
    
    # For each protein or protein-related term, count its occurrences in the file and print the count
    # The \b word boundary is used to ensure whole words are matched
    # Note: Some of the regular expressions include specific formatting or additional words for more precise matching
    count = sum(1 for match in re.finditer(r"\bhydrogenase", contents))
    print("Hidrogenases:", count)
    
    # Repeat this process for each protein or protein-related term of interest
    # ...
