#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by: Jorge Camarero Vera, Jose Jesus Gallego-Parrilla
Maintained by: Jose Jesus Gallego-Parrilla
Email: J.J.Gallego-Parrilla2@newcastle.ac.uk
License: GPL

This script uses the ProteinHtml class to extract specific protein information from a provided HTML file.
"""

# Import the ProteinHtml class from the ProteinHtml.py module
from ProteinHtml import ProteinHtml 

# This is the main entry point of the script. It checks whether this script is being run directly 
# (as opposed to being imported as a module), in which case it executes the code within this block.
if __name__ == "__main__":

    # Create an instance of the ProteinHtml class, initializing it with the name of an HTML file
    protein = ProteinHtml("chandra.html") #name of the file you want check

    # Call the getFileName method of the protein object, which returns the name of the file it was initialized with,
    # and print this file name to the console
    name = protein.getFileName()
    print(name)

    # Call the searchProtein method of the protein object, which returns the number of occurrences of a specified 
    # protein in the HTML file, and print this count to the console
    count = protein.searchProtein("Acidobacterium ailaaui")
    print(count)

    # Call the findAllProteinNames method of the protein object, which returns a list of all unique protein names found
    # in the HTML file, and print this list to the console
    found = protein.findAllProteinNames()
    print(found)
