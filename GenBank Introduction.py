# -*- coding: utf-8 -*-
"""
NCBI's databases, such as PubMed, GenBank, GEO, and many others, can be accessed via Entrez, a data retrieval system offered by NCBI. For direct access to Entrez, you can use Biopython’s Bio.Entrez module.

The Bio.Entrez.esearch() function will search any of the NCBI databases. This function takes the following arguments:

    db : The database to search. For example, this field can be nucleotide for GenBank or pubmed for PubMed.
    term: The search term for the "Query" field. You can use search tags here.

We will now demonstrate a quick search for the rbcL gene in corn (Zea mays):
"""

from Bio import Entrez
Entrez.email = "dybgreg@gmail.com"
handle = Entrez.esearch(db="pubmed", term='Exaiptasia[Organism] AND ("2000/11/01"[PDAT] : "2021/12/21"[PDAT]) ')
record = Entrez.read(handle)
record["Count"]
# Surely this value will change over time because GenBank is constantly updated

print (record["Count"])
"""
Note that when you request Entrez databases you must obey NCBI's requirements:

    For any series of more than 100 requests, access the database on the weekend or outside peak times in the US.
    Make no more than three requests every second.
    Fill in the Entrez.email field so that NCBI can contact you if there is a problem.
    Be sensible with your usage levels; if you want to download whole mammalian genomes, use NCBI's FTP.
"""
