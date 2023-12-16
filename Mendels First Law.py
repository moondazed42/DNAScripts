# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from itertools import product

k = 23  # AA  homozygous dominant
m = 20  # Aa  heterozygous
n = 20  # aa  homozygous recessive

# creates a list called population by concatenating three lists, each representing a genotype, repeated the specified number of times.
population = (['AA'] * k) + (['Aa'] * m) + (['aa'] * n)


# iterates through all possible pairs of parents in the population and calculates all possible children according to the Punnett square. 
# It uses the product function to obtain all combinations of alleles from the two parents, and then joins these alleles to form the genotype of each child.
all_children = []
for parent1 in population:
    # remove selected parent from population.
    chosen = population[:]
    chosen.remove(parent1)
    for parent2 in chosen:
        # get all possible children from 2 parents. Punnet square
        children = product(parent1, parent2)
        all_children.extend([''.join(c) for c in children])

dominants = filter(lambda c: 'A' in c, all_children)
# float for python2
print (float(len(list(dominants))) / len(all_children))
# 0.7833333
