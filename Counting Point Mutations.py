# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 01:00:53 2023

@author: Dylan

Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t)

"""

file_path = "C:\\Users\\Dylan\\Downloads\\rosalind_hamm.txt"
with open(file_path) as file:
    file_content = file.read()

seqs = file_content.split('\n')

def count_mismatches(seq1, seq2):
    # Ensure both sequences have the same length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length.")
    # Count the number of mismatches
    mismatch_count = sum(1 for base1, base2 in zip(seq1, seq2) if base1 != base2)

    return mismatch_count

def compare_sequences(sequence_list):
    # Ensure there are at least two sequences to compare
    if len(sequence_list) < 2:
        raise ValueError("At least two sequences are required for comparison.")
    # Get the length of the sequences
    sequence_length = len(sequence_list[0])
    # Iterate through the positions and count mismatches
    mismatch_counts = [count_mismatches(sequence_list[0][i], sequence_list[1][i]) for i in range(sequence_length)]
    return mismatch_counts

def tally(mutations):
    # Count the number of mismatches
    return mutations.count(1)

# Example usage
sequences = ['GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT']
mismatch_counts = compare_sequences(seqs)
count = tally(mismatch_counts)
print (mismatch_counts)
print (count)


# OR ~~~~ people who write simple code on these problems astound me. Although the upper one logs locations of mutations and that's useful.

print ([ a!=b for (a, b) in zip(seqs[0], seqs[1])].count(True))


