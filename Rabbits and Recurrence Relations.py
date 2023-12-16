# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 18:21:46 2023

In 1202, Leonardo of Pisa (commonly known as Fibonacci) considered a mathematical exercise regarding the reproduction of a population of rabbits. 

He made the following simplifying assumptions about the population:

    The population begins in the first month with a pair of newborn rabbits.
    Rabbits reach reproductive age after one month.
    In any given month, every rabbit of reproductive age mates with another rabbit of reproductive age.
    Exactly one month after two rabbits mate, they produce one male and one female rabbit.
    Rabbits never die or stop reproducing.

Fibonacci's exercise was to calculate how many pairs of rabbits would remain in one year. We can see that in the second month, the first pair 
of rabbits reach reproductive age and mate. In the third month, another pair of rabbits is born, and we have two rabbit pairs; our first pair 
of rabbits mates again. In the fourth month, another pair of rabbits is born to the original pair, while the second pair reach maturity and mate 
(with three total pairs). The dynamics of the rabbit population are illustrated in Figure 1. After a year, the rabbit population boasts 144 pairs.

Although Fibonacci's assumption of the rabbits' immortality may seem a bit farfetched, his model was not unrealistic for reproduction in a predator-free 
environment: European rabbits were introduced to Australia in the mid 19th Century, a place with no real indigenous predators for them. Within 50 years, 
the rabbits had already eradicated many plant species across the continent, leading to irreversible changes in the Australian ecosystem and turning much 
of its grasslands into eroded, practically uninhabitable parts of the modern Outback (see Figure 2). In this problem, we will use the simple idea of counting 
rabbits to introduce a new computational topic, which involves building up large solutions from smaller ones.

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the 
finite sequence (π,−2–√,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits 
from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the 
number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fnrepresents the number of rabbit pairs 
alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1to 
initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively 
larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the 
answers to smaller cases.

Given: Positive integers n≤40 and k≤5

Return: The total number of rabbit pairs that will be present after n
months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

Sample Dataset
5 3
Sample Output
19

~~~~~~~~~~~

FROM CHATGPT
This problem is a classic example of the Fibonacci sequence, but with a twist. Instead of each pair producing only one pair of offspring, 
they produce kk pairs. The recurrence relation for this problem can be expressed as follows:

F(n)=F(n−1)+k⋅F(n−2)F(n)=F(n−1)+k⋅F(n−2)

Here's a Python function to solve this problem using dynamic programming (memoization):

python

def rabbit_pairs(n, k, memo={}):
    if n == 1 or n == 2:
        return 1
    if n not in memo:
        memo[n] = rabbit_pairs(n-1, k, memo) + k * rabbit_pairs(n-2, k, memo)
    return memo[n]

# Example usage
n, k = 5, 3
result = rabbit_pairs(n, k)
print(result)  # Output: 19

In this example, n represents the number of months, and k represents the number of pairs each pair produces. 
The memo dictionary is used for memoization to avoid redundant calculations.

You can modify the values of n and k based on your input, and the function will return the total number of rabbit pairs after n months.


"""

def rabbit_pairs(n, k, den={}):
    if n == 1 or n == 2:
        return 1
    if n not in den:
        den[n] = rabbit_pairs(n-1, k, den) + k * rabbit_pairs(n-2, k, den)
    return den[n]


n, k = 6, 4
result = rabbit_pairs(n, k)
print(result)

"""

rabbit_pairs(n-1, k, memo): This part calculates the number of rabbit pairs after (n−1) months. It calls the rabbit_pairs function 
recursively with the argument n−1, representing the previous month. This part represents the number of pairs from the previous month.

k * rabbit_pairs(n-2, k, memo): This part calculates the number of pairs produced by the pairs that were reproductive (n−2) months ago. 
It multiplies k (the number of pairs each pair produces) by the number of pairs after (n−2) months. This part represents the offspring 
produced by the reproductive pairs from two months ago.

memo[n] = ... + ...: The results from the two calculations above are added together, and the total is stored in the memo dictionary at key n. 
This is a step of memoization, where we store the result for n to avoid recalculating it if it's needed again in the future.

"""