# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 21:01:35 2023

@author: Dylan
"""

fastq = "C:\\Users\\Dylan\\SandBox\\Text\\fastq.txt"
fasta = "C:\\Users\\Dylan\\SandBox\\Text\\fasta.txt"
from Bio import SeqIO
with open(fasta, 'w') as file:
    # Append content to the file
    for seq_record in SeqIO.parse(fastq, "fastq"):
        file.write('>' + seq_record.id + '\n')
        file.write(str(seq_record.seq) + '\n')