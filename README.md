# Graceful labeling + sequential labeling
Various code related to the problem of graph labelings (specifically, trees: graceful labelings and relatives of harmonious labelings)

# And introducing new beta+seq labeling!
![Tree labelings](/images/tree_labelings.png)

Notation is borrowed from [1] and [2]:
* \rho bi = \rho bigraceful
* harm = harmonious
* semt = super-edge-magic total
* seq - sequential

Beta+seq labeling - is a new labeling, which is a strict subset of both beta+ and sequential labelings (TODO: is it equal to intersection of them both?): we have 2 partitions of vertices in the tree and edges connect vertices from different partitions; to switch from beta+ to seq and vice versa we "reverse" values in one of the partitions and shift them (mod number of vertices) by appropriate value (TODO: i guess this value is easy to calculate, depending on the size of partitions). And so the conjecture is:

For every tree there exists a beta+seq labeling.

I've verified this conjecture for all trees upto and including 15 or 16 vertices (TODO).

# References:

[1]  A. Blinco, S.I. El-Zanati, C. Vanden Eynden, On the cyclic decomposition of complete graphs into almost-bipartite graphs, Discrete Math. 284 (2004) 71–81.
[2] J. A. Gallian, A dynamic survey of graph labeling, Elect. J. Combin., (2015), #DS6 - www.combinatorics.org/ojs/index.php/eljc/article/viewFile/DS6/pdf
