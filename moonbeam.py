# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 22:03:24 2024

@author: Dylan

intent: build small pipeline to analize sequences of varying length by blasting and outputing data table of summary metrics

"""

# Establishes file paths
csv_file = "C:\\filepath\\Bat Species Database Sequencing Status.csv"
output_fasta = "C:\\filepath\\moonbeam.fasta"
blast_out = "C:\\filepath\\blast_out.xml"
alignment_out = "C:\\filepath\\alignment.xml"
api = ""
email = ""

import csv
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import threading

# Function to read data from CSV and create a list of tuples | wanting to keep genus and species separate objects in case i need that later
def create_species_list(csv_file):
    data_list = []
    with open(csv_file, 'r', encoding='utf-8') as file:  # Specify encoding as utf-8
        reader = csv.reader(file)
        next(reader)  # Skip the header row if it exists
        for row in reader:
            # Assuming the second column is the genus and the third column is the species
            genus = row[1]
            species = row[2]
            data_list.append((genus, species))
    return data_list

# Function to search NCBI for accession numbers of species and create dictionary with species as key and accession list as value
def search_species(species_list):
    print ("Searching for genbank entries...")
    Entrez.email = email
    Entrez.api_key = api
    species_accessions = {}  # Dictionary to store species and their corresponding accession numbers
    total = 0  # Initialize total count
    
    for species in species_list:
        term = f'{species}[Organism]'
        handle = Entrez.esearch(db="nucleotide", term=term)
        record = Entrez.read(handle)
        count = record["Count"]
        print(f"{species}: {count}")
        accession_numbers = record["IdList"]  # List of accession numbers for the species
        if accession_numbers:
            species_accessions[f"{species}"] = accession_numbers
            total += int(count)  # Accumulate the count
        
    return species_accessions, total
    
# Function to get description from accession number        
def get_species(accession):
    Entrez.email = email
    Entrez.api_key = api
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, 'genbank')
    handle.close()
    species = record.description
    return species
    
# Function to identify if each species has full genomes        
def search_genomes(species_list):
    print("Searching for whole alignments...")
    Entrez.email = email
    Entrez.api_key = api
    for species in species_list:
        handle = Entrez.esearch(db="genome", term=f'{species}[Organism]')
        record = Entrez.read(handle)
        count = record["Count"]
        print(f"{species}: {count}")
        
# Function to fetch fastas from accession numbers
def fetch_fasta_from_ncbi(accession):
    try:
        Entrez.email = email
        Entrez.api_key = api
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
        if len(record.seq) >= 2000:  # Check sequence lengths, filter out any less than 1000bp
            return record
        else:
            return None
    except Exception:
        return None

# Function to create a lookup table of compiled entires, as well as compile a local FASTA file for all accession numbers
def create_fasta_lookup_table(species_accessions, output_fasta):
    fasta_accession_sequences = {}  # Dictionary to store accession numbers and their corresponding FASTA sequences
    print ("Creating FASTA lookup table...")
    with open(output_fasta, "w") as out_handle:
        for species, accession_numbers in species_accessions.items():
            fasta_sequences = []
            for accession in accession_numbers:
                fasta_record = fetch_fasta_from_ncbi(accession)
                if fasta_record:
                    fasta_sequences.append(fasta_record)
                    SeqIO.write(fasta_record, out_handle, "fasta")
                    fasta_accession_sequences[accession] = fasta_record.seq
                else:
                    continue
                print (species, accession, fasta_record)
    return fasta_accession_sequences
    
# Functon to blast sequence and pull alignment from blast results
def blast_sequence(seq_object, accession, db, blast_out, alignment_out, timeout_seconds, species):
    Entrez.email = email
    Entrez.api_key = api
    try:
        # Define a function to perform the blast search
        def perform_blast():
            result_handle = NCBIWWW.qblast('blastn', db, seq_object, megablast=True, format_type="XML")
            response = result_handle.read()
            result_handle.close()
            with open(blast_out, "w") as out_handle:
                out_handle.write(response)
            print(f"BLAST results appended to {blast_out}")

        # Create a thread for the blast search
        blast_thread = threading.Thread(target=perform_blast)
        blast_thread.start()
        
        # Wait for the thread to finish or timeout
        blast_thread.join(timeout_seconds)
        
        # Process the BLAST results
        E_VALUE_THRESH = 1.00
        with open(blast_out) as result_handle, open(alignment_out, "a") as alignment_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            alignment_handle.write(f"Hit found for {species}, sequence accession ID: {accession}\n")
                            alignment_handle.write("****Alignment****\n")
                            alignment_handle.write(f"Sequence: {alignment.title}\n")
                            alignment_handle.write(f"Length: {alignment.length}\n")
                            alignment_handle.write(f"E Value: {hsp.expect}\n")
                            alignment_handle.write(f"Query: {hsp.query[:75]}...\n")
                            alignment_handle.write(f"Match: {hsp.match[:75]}...\n")
                            alignment_handle.write(f"Subject: {hsp.sbjct[:75]}...\n")
                            alignment_handle.write("\n")  # Separator between alignments
                            print(f"Hit found for {species}, sequence accession ID: {accession}")
                            print("****Alignment****\n")
                            print(f"Sequence: {alignment.title}\n")
                            print(f"Length: {alignment.length}\n")
                            print(f"E Value: {hsp.expect}\n")
                            print(f"Query: {hsp.query[:75]}...\n")
                            print(f"Match: {hsp.match[:75]}...\n")
                            print(f"Subject: {hsp.sbjct[:75]}...\n")
                            print("\n")  # Separator between alignments
            print(f"Alignments written to {alignment_out}")
            
            # Check if the thread has timed out
            if blast_thread.is_alive():
                print(f"BLAST search for sequence {accession} took longer than {timeout_seconds} seconds and timed out. Moving on to the next sequence.")
                blast_thread.cancel()
                return  # Exit the function if the thread timed out

    except Exception as e:
        print(f"An error occurred during BLAST search for {accession}: {e}")

# Prepares sequencing list from csv file of species lists supplied by Bat1K
species_list = [" ".join(t) for t in create_species_list(csv_file)]
genomes = search_genomes(species_list)
species_accessions, total = search_species(species_list)  # Unpack the tuple returned by search_species_in_genbank
print(f"Total FASTA entries: {total}")  # Print the total count of FASTA entries

# Creates lookup table and local FASTA file
fasta_accession_sequences = create_fasta_lookup_table(species_accessions, output_fasta)
print (f"Total Sequences collected: {len(fasta_accession_sequences)}")

# iterates through sequences in the lookup table and blasts sequences against virus database, skips sequence if blast takes longer than n seconds
for accession, seq in fasta_accession_sequences.items():
    species = get_species(accession)
    print(f"Blasting {species}, accessionID {accession} with sequence:")
    print(seq)
    db = 'nt_viruses'
    timeout_seconds = 600  # Timeout set to 300 seconds (5 minutes)
    blast_sequence(seq, accession, db, blast_out, alignment_out, timeout_seconds, species)
