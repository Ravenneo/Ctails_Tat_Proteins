#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 1 12:12:05 2020

This module contains the ProteinHtml class, which is used to analyze protein information in an HTML file.
"""

# Import the re module for regular expression operations
import re


class ProteinHtml:
  # Private class variables to hold the name of the file and its content
  __fileName = None
  __contentFile = None

  # Constructor method that initializes a new instance of the class
  def __init__(self, fileName):
    self.__fileName = fileName  # Assign the provided file name to the __fileName variable
    with open(fileName) as f:  # Open the file and assign its content to the __contentFile variable
      self.__contentFile = f.read()

  # Method to get the name of the file
  def getFileName(self):
    return self.__fileName

  # Method to count the number of occurrences of a specific protein in the file
  def searchProtein(self, proteinName):
    # Construct a regular expression that matches the specified protein in the HTML structure
    regexBegin = "<tr style = \"background:#[\w\d]+\"><td>WP_[\d]+\.1<\/td><td colspan = 5>(.*"
    regexEnd = ".*)<\/td><\/tr>"
    regex = regexBegin + proteinName + regexEnd

    # Count the number of matches of the regex in the content of the file
    count = sum(1 for match in re.finditer(r"{}".format(regex),
                                           self.__contentFile))
    return count

  # Method to find all unique protein names in the file
  def findAllProteinNames(self):
    allProteins = {}  # Dictionary to hold each protein name and its count of occurrences

    # Construct a regular expression that matches protein names in the HTML structure
    regex = "<tr style = \"background:#[\w\d]+\"><td>WP_[\d]+\.1<\/td><td colspan = 5>.*\[(.*)\]<\/td><\/tr>"

    # Iterate over all matches of the regex in the content of the file
    for match in re.finditer(r"{}".format(regex), self.__contentFile):
      protein = match.group(1)  # Extract the protein name from the match

      # If the protein is already in the dictionary, increment its count;
      # otherwise, add it to the dictionary with a count of 1
      if (protein in allProteins):
        allProteins[protein] += 1
      else:
        allProteins.update({protein: 1})

    return allProteins
