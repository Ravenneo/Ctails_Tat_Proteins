#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script allows the user to interactively choose between two options: displaying protein ID in the console or exporting protein IDs to a CSV file.

Created on Mon Jul 13 19:30:19 2020
__author__ = 'Jose Jesus Gallego-Parrilla'
__license__ = "GPL"
__maintainer__ = "Jose Jesus Gallego-Parrilla"
__email__ = "J.J.Gallego-Parrilla2@newcastle.ac.uk"
"""

from ProteinHtmlID import ProteinHtmlID
import pandas as pd

# Interactively ask the user to choose an option
while True:
    try:
        answer = int(input("Press 1 to see protein ID in console \nPress 2 to export protein CSV list \nChoose="))
        if answer == 1 or answer == 2:
            break
    except ValueError:
        pass
    print("Sorry, not what I was expecting \nTry again")

if answer== 1:
    # If the user chooses option 1, display the protein ID in the console
    protein = ProteinHtmlID("chandra.html") #name of the file you want check
    name = protein.getFileName()
    print(name)
    
    count = protein.searchProtein("Acidobacterium ailaaui")
    print(count)
    
    found = protein.findAllProteinNames()
    print(found)
elif answer== 2:
    # If the user chooses option 2, export the protein IDs to a CSV file
    wp_num = ProteinHtmlID("chandra.html")
    found = wp_num.findAllProteinNames()
    
    wp_num = []
    for elem in found:
        wp_num.append(elem)
      
    data = {'ID' : wp_num}
    dataframe = pd.DataFrame(data)
    input_console = input("Give the file a name that ends in .csv:\n")
    dataframe.to_csv(input_console, index=False) #name for exported document
    print(dataframe)
