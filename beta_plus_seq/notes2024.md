TODO: retrieve notes from cycle-double-cover project


TODO: ? bijection between beta and seq labelings
    TODO: number of sequential labelings for all graphs (not only trees) with n vertices is equal to n!
        TODO: prove using Wilf-Zeilberger method
            or just prove it recursively
    TODO: is it actually the same set of graphs?
        for n=1,2,3 this is true


TODO: is beta+seq equal to intersection of them both?
    TODO: bipartite graphs
    TODO: complete bipartite graphs
    TODO: crown graphs, formed from complete bipartite graphs by removing the edges of a perfect matching
    TODO: even cycle graphs
    TODO: Every planar graph whose faces all have even length is bipartite
    TODO: bipartite hamiltonian cycle graphs
    TODO: Hypercube graphs, partial cubes, and median graphs
        In these graphs, the vertices may be labeled by bitvectors, in such a way that two vertices are adjacent if and only if the corresponding bitvectors differ in a single position

    β+seq labeling - is a new labeling, which is a strict subset of both β+ and sequential labelings (?)
    we have 2 partitions of vertices in the tree and edges connect vertices from different partitions
    to switch from β+ to seq and vice versa we "reverse" values in one of the partitions and shift them (mod number of vertices) by appropriate value


TODO: calculate shift value between beta+ and seq
    ? depending on the size of partitions?


TODO: S_{n, 2} family of trees
    S3:
    T = {0->1, 1->2, 0->3, 3->4, 0->5, 5->6}
    β+ labels: 0 6 2 5 4 3 1
    seq labels: 2 6 0 5 4 3 1
    β+ edges: 0->6, 6->2, 0->5, 5->4, 0->3, 3->1 (which means that we have every difference between 6 and 1)
    seq edges: 2->6, 6->0, 2->5, 5->4, 2->3, 3->1 (which means that we have every sum between 4 and 9)

    S4:
    ...

    TODO: there was some magic (magical feeling to the labelings) also, right?
