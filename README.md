# Graceful labeling + sequential labeling + super-edge-graceful labeling (todo)

(btw, WARNING: LICENSE is not applicable to any of freetree.h files; consult [http://pallini.di.uniroma1.it](http://pallini.di.uniroma1.it) for licensing)

[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)


Various code related to the problem of graph labelings (specifically, trees: graceful labelings and relatives of harmonious labelings)

# And introducing new β+seq labeling!
![Tree labelings](/diagram/tree_labelings.png)

Notation is borrowed from [1] and [2]:
* ρ++ bi = ρ++ bigraceful
* black dotted border means that this labeling has known families of counterexamples and it is conjectured that there are no other counterexamples
* red colour means that this labeling doesn't work for all trees and we don't have a full characterization of counterexamples
* yellow colour means that it's a theorem
* bold typeface means that the labeling is important or famous

β+seq labeling - is a new labeling, which is a strict subset of both β+ and sequential labelings (TODO: is it equal to intersection of them both?): we have 2 partitions of vertices in the tree and edges connect vertices from different partitions; to switch from β+ to seq and vice versa we "reverse" values in one of the partitions and shift them (mod number of vertices) by appropriate value (TODO: i guess this value is easy to calculate, depending on the size of partitions). And so the conjecture is:

For every tree there exists a β+seq labeling.

I've verified this conjecture for all trees upto and including 15 or 16 vertices (TODO).

It's also very easy to see that every alpha labeling is also a β+seq labeling (by the same argument: reverse the values in one of the halves; and the shift is 0)

Example:
* T = {0->1, 1->2, 0->3, 3->4, 0->5, 5->6}
* β+ labels: 0 6 2 5 4 3 1
* seq labels: 2 6 0 5 4 3 1
* β+ edges: 0->6, 6->2, 0->5, 5->4, 0->3, 3->1 (which means that we have every difference between 6 and 1)
* seq edges: 2->6, 6->0, 2->5, 5->4, 2->3, 3->1 (which means that we have every sum between 4 and 9)

# Other new results (or alternatively, it's also interesting to note the following:)
* There is also a very interesting any-β+ labeling (when we fix an arbitrary set of values for edges and still ask for properties of β+ that values of vertices in one of the partitions are locally more than values of another partition). It looks like there is only 1 affine family of counterexamples: T = {0->1, 1->2, 0->3, 3->4, 0->5, 5->6}; E = {1, 2, 3, 4, 6, 7} * c, c >= 1 (and also a question: does there exist an any-sequential labeling? Maybe also an any-β+seq labeling?)
* already mentioned any-β+ labeling
* ρ++ bigraceful labeling is the same as σ++ labeling (you just need to reverse values in one the partitions as in β+seq)
* and actually all other new arrows and equivalences in the diagram (e. g., odd harmonious == odd-graceful, harmonious -> ρ bigraceful & c.)
* number of graceful labelings for trees of fixed number of vertices is always divisible by 4 (TODO: we get a factor of 2 by reversing all labels and another factor by 2 by applying a proof by moving the edge with label 1) (oh, actually this was also recently independently proven in [8])
* number of sequential labelings for all graphs (not only trees) with n vertices is equal to (n-1)! (TODO: proof using Wilf-Zeilberger method); which is actually the same number of graceful labelings for graphs with n vertices. Maybe there do exist some interesting bijection between graceful and sequential labelings?
* and maybe, but no intuition or proof, with n not equal to 3, the number of alpha-labelings for all graphs is always an odd number (verified for n <= 14)
* suppose we study a question of existence of alpha labeling of a tree. Sometimes a tree doesn't have it because of the parity condition (TODO). Can this obstruction also be a problem for β+ labeling? It turns out that no, it can't (TODO).
* there's a hypothetical characterization of all treees which are not 0-ubiquitous [3]. We know that maximally valued edge is always connected to the 0 vertex. Can we characterize the forbidden places of this edge in trees? It looks like it's possible, and in addition to 2 interesting families in [3] of forbidden placings of this edge, we get only 4 (or 5?) new families (TODO). It also seems that we do get something similar for sequential labeling, although the number of families of counterexamples looks huge or infinite at the moment (TODO).
* conjecturing here that elegant labeling exists for all trees, except for 1 family of counterexamples (TODO)
* the inductive constructions by Koh, Rogers and Tan [4, 5, and also see 6, 7] and the transfer method both are useful for determining the possible placings of the maximally valued edge (in case of so called delta construction we have 2 trees, A and B, and possible placings of maximal-value edge in both of them and from this information we can derive the possible placings of maximal-value edge in their delta and delta+ products), but i haven't tried to code this yet.
* Code for graceful tree packing conjecture, where every tree recieves a similar-to-graceful labeling (TODO). And more!
# References:

[1]  A. Blinco, S.I. El-Zanati, C. Vanden Eynden, On the cyclic decomposition of complete graphs into almost-bipartite graphs, Discrete Math. 284 (2004) 71–81.

[2] J. A. Gallian, A dynamic survey of graph labeling, Elect. J. Combin., (2015), #DS6 - www.combinatorics.org/ojs/index.php/eljc/article/viewFile/DS6/pdf

[3] F. Van Bussel, 0-Centred and 0-ubiquitously graceful trees, Discrete Math. 277: 193 − 218 (2004).

[4] K. M. Koh, D. G. Rogers and T. Tan, “On Graceful Trees”, Nanta Mathematica, Vol. 10, No. 2, pp. 27–31, 1977.

[5] K. M. Koh, D. G. Rogers and T. Tan, “Products of Graceful Trees”, Discrete Mathematics, Vol. 31, pp. 279–292, 1980.

[6] M. Mavronicolas, L. Michael, "A substitution theorem for graceful trees and its applications", Discrete Mathematics, Vol. 309, pp. 3757–3766, 2009.

[7] M. C. Superdock, The graceful tree conjecture: a class of graceful diameter-6 trees, arXiv preprint arXiv:1403.1564, (2014).

[8] D. Anick, "Counting graceful labelings of trees: A theoretical and empirical study", Discrete Applied Mathematics, Vol. 198, pp. 65-81, 2016.
