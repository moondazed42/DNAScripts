# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:44:30 2023

@author: Dylan

Counting DNA Nucleotides

"""
print ('imput sequence:')
pw = input('')

dna = {}
for base in pw:
    if base in dna:
        dna[base] += 1
    else:
        dna[base] = 1
for base, count in dna.items():
    print (base, count)
    
print (dna)