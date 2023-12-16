# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 21:01:35 2023

@author: Dylan
"""

fastq = "C:\\Users\\Dylan\\SandBox\\Text\\fastq.txt"
fasta = "C:\\Users\\Dylan\\SandBox\\Text\\fasta.txt"

from Bio import SeqIO       
        
with open(fastq) as input_data, open(fasta, 'w') as output_data:
    SeqIO.convert(input_data, 'fastq', output_data, 'fasta' )