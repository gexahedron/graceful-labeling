β+seq labeling - is a new labeling, which is a strict subset of both β+ and sequential labelings
we have 2 partitions of vertices in the tree and edges connect vertices from different partitions
to switch from β+ to seq and vice versa we "reverse" values in one of the partitions and shift them (mod some number) by appropriate value
so formula is (shift1 - label) mod shift2


TODO: retrieve notes from cycle-double-cover project


TODO: ? bijection between beta and seq labelings
    TODO: number of sequential labelings for all graphs (not only trees) with n vertices is equal to (n-1)!
        TODO: prove using Wilf-Zeilberger method
            or just prove it recursively

    DONE: is it actually the same set of graphs?
        for n=2,3,4 this is true
        no, breaks on n=5
    DONE: is it just same graphs, without counting them?
        no, 6 vertices, cycle of length 5
        it's not graceful, but sequential
        edges {(2, 4), (0, 3), (1, 4), (0, 2), (1, 3)}
        sums: 2,3,4,5,6
    TODO/DONE: is graceful graphs a subset of sequential graphs?
        no, here are counts of isomorphic graphs, (n, beta vs seq):
        2 1 1
        3 1 1
        4 3 3
        5 5 5
        6 12 14
            seq [(0, 3), (0, 2), (1, 4), (1, 3), (2, 4)]
            seq [(0, 5), (1, 3), (1, 2), (2, 4), (3, 4)]
        7 37 34 (so we even have more beta graphs than seq)
            beta [(0, 6), (0, 5), (1, 5), (1, 4), (2, 3), (4, 6)]
            beta [(0, 6), (0, 5), (1, 4), (1, 3), (2, 6), (3, 4)]
            beta [(0, 6), (0, 5), (1, 3), (2, 6), (2, 5), (3, 4)]
        8 112 109
            beta [(0, 7), (0, 6), (1, 6), (1, 3), (2, 5), (3, 7), (4, 5)]
                01367 245
            beta [(0, 7), (0, 6), (1, 6), (2, 5), (2, 4), (3, 7), (4, 5)]
                01367 245
            beta [(0, 7), (0, 6), (1, 4), (1, 3), (2, 7), (2, 6), (3, 4)]
                0267 134
        9 340 337
            beta [(0, 8), (0, 7), (0, 6), (0, 3), (1, 5), (2, 4), (3, 8), (6, 7)]
                03678 15 24
            beta [(0, 8), (0, 7), (1, 7), (1, 5), (2, 5), (2, 3), (3, 8), (4, 6)]
                0123578 46
            beta [(0, 8), (0, 7), (0, 3), (1, 5), (2, 8), (2, 7), (2, 3), (4, 6)]
                02378 15 46
            beta [(0, 8), (0, 7), (0, 2), (1, 5), (1, 4), (2, 8), (2, 7), (4, 5)]
                0278 145
            beta [(0, 8), (0, 7), (0, 2), (1, 5), (2, 8), (2, 7), (3, 6), (7, 8)]
                0278 15 36
            beta [(0, 8), (0, 7), (1, 5), (2, 8), (2, 7), (3, 6), (3, 4), (4, 6)]
                0278 15 346
            seq [(0, 8), (0, 7), (0, 5), (1, 3), (1, 2), (2, 4), (3, 6), (4, 6)]
                0578 12346
            seq [(0, 7), (0, 5), (1, 8), (1, 3), (2, 6), (2, 4), (3, 8), (4, 6)]
                057 138 246
            seq [(0, 7), (0, 4), (1, 5), (2, 8), (2, 6), (2, 3), (3, 6), (4, 7)]
                047 15 2368
        so, the difference between the sets is not that big somehow
        and it's mostly (1 exception) about when we have more that 1 connected component

    TODO/DONE: but it's interesting, if we compare the degree sequences only
        then it's true that beta are subset of seq, at least upto 12 vertices including
        they are even exactly same for 1,2,3,4,5,7,8,9,11,12 vertices
        so only strict subset for 6 and 10


DONE: is beta+seq equal to intersection of them both?
    no,
        seq [(0, 5), (1, 3), (1, 2), (2, 4), (3, 4)]
        there's no beta-labeling of this graph
        it's disconnected though
        but, it's bipartite and planar; 4-cycle + edge
TODO: is beta+seq equal to intersection of them both, on CONNECTED graphs?
    let's focus on bipartite connected graphs, where one of beta+ or seq works
    TODO: all bipartite graphs
    TODO: complete bipartite graphs
    TODO: crown graphs, formed from complete bipartite graphs by removing the edges of a perfect matching
    TODO: K_{m,n} graphs, where one of beta+ or seq works
    TODO: even cycle graphs
        NOTE: C6: no beta+ labeling - is it a parity problem?
    TODO: Every planar graph whose faces all have even length is bipartite
    TODO: bipartite hamiltonian cycle graphs
    TODO: Hypercube graphs, partial cubes, and median graphs
        In these graphs, the vertices may be labeled by bitvectors, in such a way that two vertices are adjacent if and only if the corresponding bitvectors differ in a single position

    CONJ: all simple connected bipartite graphs have beta+ labeling, except for 1 family:
        when n % 4 == 2, and gcd(degrees) = 2
            TODO: probably because of parity problems
        works at least until 12 edges including


TODO: calculate shift value between beta+ and seq
    maybe only depending on the size of partitions?


TODO: S_{n, 2} family of trees
    S3:
    T = {0->1, 1->2, 0->3, 3->4, 0->5, 5->6}
    β+ labels: 0 6 2 5 4 3 1
    seq labels: 2 6 0 5 4 3 1
    β+ edges: 0->6, 6->2, 0->5, 5->4, 0->3, 3->1 (which means that we have every difference between 6 and 1)
    seq edges: 2->6, 6->0, 2->5, 5->4, 2->3, 3->1 (which means that we have every sum between 4 and 9)

    S4:
    ...

    TODO: there was some "magic" (magical feeling to the labelings) also, right?
        i think this was more about beta+ labelings
        because beta+seq are a bit weirder, they seem to have period 3 pattern
        but there's a lot of regularity in the solutions anyway, which is interesting
        maybe there's some kind of recursion going on

    3: [6 5 3] + [2 4 1]
    6: [12 11 5 7 9 8] + [2 10 1 4 3 6]
    9: 4 options, but unordered they are same
        [18 17 7 10 9 13 11 15 14] [2 16 1 6 4 5 8 3 12]
        or
        [18 17 13 15 7 9 11 14 10] [2 16 1 12 3 4 5 6 8]
        or
        [18 7 9 10 14 11 15 13 17] [16 2 5 4 6 8 3 12 1]
        or
        [18 14 15 7 10 11 13 9 17] [16 2 12 3 4 6 5 8 1]

        ? which one to choose? there are 6 solutions
        let's look at the cycle decompositions
        [18, 17, 15, 14, 13, 11, 10, 9, 7] [2, 16, 3, 12, 5, 8, 6, 4, 1]
            [18, 17, 7, 10, 9, 13, 11, 15, 14] [2, 16, 1, 6, 4, 5, 8, 3, 12]
        [18, 17, 15, 14, 13, 11, 10, 9, 7] [2, 16, 12, 6, 1, 5, 8, 4, 3]
            [18, 17, 13, 15, 7, 9, 11, 14, 10] [2, 16, 1, 12, 3, 4, 5, 6, 8]
        [18, 17, 15, 14, 13, 11, 10, 9, 7] [16, 1, 3, 6, 12, 8, 4, 5, 2]
            [18, 7, 9, 10, 14, 11, 15, 13, 17] [16, 2, 5, 4, 6, 8, 3, 12, 1]
        [18, 17, 15, 14, 13, 11, 10, 9, 7] [16, 1, 12, 2, 5, 6, 4, 8, 3]
            [18, 14, 15, 7, 10, 11, 13, 9, 17] [16, 2, 12, 3, 4, 6, 5, 8, 1]

        [18, 17, 16, 14, 13, 11, 10, 8, 7] [15, 2, 4, 5, 12, 9, 6, 3, 1]
            [18, 8, 14, 11, 17] [15, 3, 5, 9, 2]
            [16, 13, 7, 10] [4, 12, 1, 6]
        [18, 17, 16, 14, 13, 11, 10, 8, 7] [15, 5, 1, 12, 4, 6, 9, 2, 3]
            [18, 7, 13, 10, 16] [15, 3, 4, 9, 1]
            [17, 14, 8, 11] [5, 12, 2, 6]

    unordered:
    3: 3,5,6 + 1,2,4
        0..o.oo
    6: 5,7,8,9,11,12 + 1,2,3,4,6,10
        0....o.ooo.oo
    9: 7,9,10,11,13,14,15,17,18 + 1,2,3,4,5,6,8,12,16
        0......o.ooo.ooo.oo
    looks like there is period 4 pattern in the values of numbers in the tree
    DONE: test 12
        looks like this approach breaks on 12

    S_{12, 2}

        False [0, 24, 23, 22, 21, 20, 17, 16, 15, 12, 11, 10, 9, 6, 4, 8, 18, 19, 13, 14, 2, 7, 5, 3, 1] [8, 24, 23, 22, 21, 20, 17, 16, 15, 12, 11, 10, 9, 2, 4, 0, 14, 13, 19, 18, 6, 1, 3, 5, 7] 24 8
        ...
        False [0, 24, 23, 21, 20, 19, 18, 16, 14, 13, 12, 11, 9, 22, 1, 17, 5, 2, 15, 6, 7, 8, 4, 10, 3] [8, 24, 23, 21, 20, 19, 18, 16, 14, 13, 12, 11, 9, 10, 7, 15, 3, 6, 17, 2, 1, 0, 4, 22, 5] 24 8
        ...

        0
        24, 23, 22, 21, 20, 17, 16, 15, 12, 11, 10, 9,
         6,  4,  8, 18, 19, 13, 14,  2,  7,  5,  3, 1

        0........oooo..ooo..ooooo

        0
        24, 23, 21, 20, 19, 18, 16, 14, 13, 12, 11, 9,
        22,  1, 17,  5,  2, 15,  6,  7,  8,  4, 10, 3
        0........o.oooo.o.oooo.oo
