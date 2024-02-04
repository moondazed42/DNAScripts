# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 22:00:20 2023

@author: Dylan
"""

file_path = ""

with open(file_path) as file:
    file_content = file.read()
seqs = file_content.split('\n')
s = seqs[0]
t = seqs[1]
index = 0

# List to store the indices of occurrences
repeats = []

# Loop to find all occurrences
while True:
    # Find the next occurrence of the smaller string
    index = s.find(t, index)
    # If no more occurrences are found, break out of the loop
    if index == -1:
        break
    # Append the index to the list
    repeats.append(index+1)
    # Move the starting index for the next search to avoid infinite loop
    index += 1
# Display the list of occurrences
for i in repeats:
    print(i)