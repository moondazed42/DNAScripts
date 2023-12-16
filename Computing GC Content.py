# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 22:03:35 2023

@author: Dylan
"""

# Specify the path to your text file
file_path = "C:\\Users\\Dylan\\Downloads\\rosalind_gc.txt"
# Open the file in read mode

with open(file_path) as file:
    # Read the contents of the file
    file_content = file.read()

# Initialize dictionaries to store sequences and GC content
sequence_dict = {}
GC_dict = {}

# Split the file content into sequences
sequences = file_content.split('>')

# Iterate through sequences
for seq in sequences:
    
    # Skip empty sequences
    if not seq:
        continue
    
    # Split each sequence into lines
    lines = seq.strip().split('\n')
    # The first line is the header, and the rest are the DNA sequence
    header = lines[0]
    dna_sequence = ''.join(lines[1:])
    # Add the entry to the dictionary
    sequence_dict[header] = dna_sequence

# Print the resulting dictionary
for key, value in sequence_dict.items():
    A = T = C = G = 0 
    for base in value:
        if base == 'A':
            A += 1
        elif base == 'T':
            T += 1
        elif base == "C":
            C += 1
        elif base == "G":
            G += 1

    # Calculate GC content
    total_bases = A + T + C + G
    GC_content = (C + G) / total_bases if total_bases > 0 else 0

    # Add the GC content to the dictionary
    GC_dict[key] = GC_content

max_key = max(GC_dict, key=GC_dict.get)

# Print the result
print(f'{max_key}\n{GC_dict[max_key] * 100:.7f}')
