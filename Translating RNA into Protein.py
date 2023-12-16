# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 15:06:24 2023

@author: dgregory
"""

import python_codon_tables as pct
table = pct.get_codons_table("9606", replace_U_by_T=False)
print ('Enter sequence to translate:')
sequence = input()
codons = {}

for i in table:
    codons[i] = []
    for n in table[i]:
        codons[i].append(n)

aminoseq = ""
codon_list = [sequence[i:i+3] for i in range(0,len(sequence),3)]

for triplet in codon_list:
    for amino, codon_values in codons.items():
        if triplet in codon_values:
            aminoseq += amino
print ("Amino Acid Sequence:", aminoseq)