// k-beta-seq labelings; palindrome schemes
// labeling conventions:
// |V| = n, |E| = n - 1; m >= n
// E \in [1, ..., m-1], m-1 >= |E|, m-1 = max{E}
// V \in [0, 1, ..., m-1], min{V} = 0, max{V} = m-1
// m < maxn
// Generalising beta-plus - two partitions, locally one of them always bigger
// Labeling works for caterpillars quite straightforwardly, so we just skip them

// RESULTS
// T = {0->1, 1->2, 0->3, 3->4, 0->5, 5->6}
// E = {1, 2, 3, 4, 6, 7} * k, k is a natural number
// verified upto n=15, depth=8

#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include "freetree.h"

using namespace std;

// TODO: use conventions for loop variables
// v, u - for vertices; e - for edges
// TODO: needs refactoring; put everything into classes and structs and namespaces;
// Or maybe just split everything into blocks commented with their functions-experiments;
// TODO: some of these ints are actually size_t
const size_t maxn = 30;
size_t n, m;
vector<size_t> gr[maxn];

int vertexLabels[maxn];
int edgeLabels[maxn];
vector<vector<int>> betaLabels;
vector<vector<int>> betaSchemes;
vector<int> curBetaLabel;
vector<int> curBetaScheme;

int parent[maxn];
vector<int> children[maxn];

int type[maxn];
int deg[maxn];

int link[maxn];
bool isCopy[maxn];
int revLinkCount[maxn];

bool usedVertex[maxn];
bool usedEdge[maxn];

int diam;
int branch[maxn];
int distTo0[maxn];
int diams[maxn];

bool hasAlpha;

size_t markedVertex[2]; //rename!
bool res;

int kol;

clock_t beginClock, curClock, prevClock, endClock;

//k-beta+ specific structures
int needEdge[maxn];
const int maxdepth = 0;
size_t excluded[maxdepth + 1];
vector<set<size_t>> edgeSets;
//set<size_t> edgeSet;

void countDiam() {
	for (size_t i = 1; i < n; ++i) {
		if (parent[i] == 0) {
			branch[i] = i;
			diams[i] = 1;
		} else {
			branch[i] = branch[parent[i]];
        }
    }
	distTo0[0] = 0;
	for (size_t i = 1; i < n; ++i)
		distTo0[i] = distTo0[parent[i]] + 1;
	for (size_t i = 1; i < n; ++i)
		diams[branch[i]] = max(diams[branch[i]], distTo0[i]);
	diam = diams[gr[0][0]] + diams[gr[0][1]];
}

bool checkTrees(int v1, int v2) {
	if (deg[v1] != deg[v2])
		return false;
	for (size_t i = 1; i < gr[v1].size(); ++i) {
		if (!(checkTrees(gr[v1][i], gr[v2][i]))) {
			return false;
        }
    }
	return true;
}

void markAsCopy(int v1) {
	isCopy[v1] = true;
	for (size_t i = 1; i < gr[v1].size(); ++i) {
		markAsCopy(gr[v1][i]);
    }
}

void buildLinks() {
	for (size_t i = 0; i < n; ++i) {
		link[i] = -1;
		isCopy[i] = false;
	}
	for (size_t center = 0; center < n; ++center) {
		for (size_t ii = 0; ii < gr[center].size() - 1; ++ii) {
			size_t neib1 = gr[center][ii];
			if (neib1 < center)
				continue;
			int neib2 = gr[center][ii + 1];
			if (checkTrees(neib1, neib2)) {
				link[neib2] = neib1;
				if (!isCopy[neib2]) {
					markAsCopy(neib2);
                }
			}
		}
	}

	for (size_t i = 0; i < n; ++i)
		revLinkCount[i] = 0;
	for (size_t i = n - 1; i--;) { // warning: specific decreasing for with size_t
		if (link[i] != -1) {
			revLinkCount[link[i]] = revLinkCount[i] + 1;
        }
    }
}

bool gen(size_t v) {
	if (v == n) {
        cout << endl;
        for (size_t i = 0; i < n; ++i) {
            cout << vertexLabels[i] << " ";
        }
        cout << endl;
        ++kol;
        return false;//true;
	} else {
        bool isMarked = ((v == markedVertex[0]) || (v == markedVertex[1]));
		int lowVal = 1;
		int hiVal = m - 2;
		int startVal = lowVal;
		int endVal = hiVal;
		int dir = 1;

		if (!isMarked) { // k-beta+
        //if (false) {     // k-beta 
			if (type[v] == 0) {
				if (parent[v] != -1)
					hiVal = vertexLabels[parent[v]] - 1;
				startVal = lowVal;
				endVal = hiVal;
				dir = 1;
				if (startVal > endVal) {
					return false;
                }
			} else {
				if (parent[v] != -1)
					lowVal = vertexLabels[parent[v]] + 1;
				startVal = hiVal;
				endVal = lowVal;
				dir = -1;
				if (startVal < endVal) {
					return false;
                }
			}
		}

		for (int i = startVal; i != endVal + dir; i += dir) {
			if (usedVertex[i] && !isMarked) {
				continue;
			} else {
				if (!isMarked) {
					vertexLabels[v] = i;
                }

				if (parent[v] != -1) {
					edgeLabels[v] = abs(vertexLabels[v] - vertexLabels[parent[v]]);
					if (!needEdge[edgeLabels[v]] || usedEdge[edgeLabels[v]]) {
						if (isMarked) {
							break;
                        } else {
							continue;
                        }
                    }

					if (link[v] != -1) {
						if (edgeLabels[v] > edgeLabels[link[v]]) {
							continue;
                        }
                    }
					if (edgeLabels[v] < revLinkCount[v]) {
						continue;
                    }
                }
					
                usedEdge[edgeLabels[v]] = true;
                usedVertex[vertexLabels[v]] = true;

    			if (gen(v + 1))
	    			return true;

				usedEdge[edgeLabels[v]] = false;
				usedVertex[vertexLabels[v]] = false;

				if (isMarked) {
					break;
                }
			}
		}
	}
	return false;
}

bool checkCatterpillarness() {
    for (size_t v = 0; v < n; ++v) {
        if (deg[v] > 2) {
            int kolDegNot1 = 0;
            for (size_t i = 0; i < gr[v].size(); ++i) {
                if (deg[gr[v][i]] > 1) {
                    ++kolDegNot1;
                }
            }
            if (kolDegNot1 > 2) {
                return false;
            }
        }
    }
	return true;
}

void genEdgeSets(int depth) {
    if (depth == -1) {
       set<size_t> s;
       for (size_t i = 1; i < excluded[maxdepth]; ++i)
           s.insert(i);
       for (size_t i = 0; i < maxdepth; ++i)
           s.erase(excluded[i]);
       edgeSets.push_back(s);
       return;
    }
    for (size_t i = excluded[depth + 1] - 1; i >= 1; --i) {
        excluded[depth] = i;
        genEdgeSets(depth - 1);
    }
}

int main(int /*argc*/, char** /*argv*/)
{
	//bool tempCode = freopen("test.txt", "w", stdout);
	beginClock = clock();
	prevClock = beginClock;
	srand (time(NULL));

    static const bool USE_CONCRETE_TREE = false;
    if (USE_CONCRETE_TREE) {
        /*totalAmount = 1;
        //const int V = 7;
        //int curTree[V + 1] = {V, -1, 0, 1, 0, 3, 0, 5};
        const int V = 10;
        int curTree[V + 1] = {V, -1, 0, 1, 2, 0, 4, 5, 0, 7, 8};
        for (size_t i = 0; i < V + 1; ++i)
            treeDatabase[0][i] = curTree[i];
        for (size_t i = 0; i < V - 1; ++i) {
            size_t val;
            cin >> val;
            edgeSet.insert(val);
        }*/
    } else {
        size_t razmBegin = 8; // rename
        size_t razmEnd = 9; // rename
    	for (size_t i = razmBegin; i <= razmEnd; ++i)
	    	freeTree(i);
    	cerr << "#trees=" << totalAmount << endl;
    }

	size_t prevn = 0; // rename slightly
	for (size_t ttt = 0; ttt < totalAmount; ++ttt) { // rename
        size_t tt = ttt; // rename
        n = static_cast<size_t>(treeDatabase[tt][0]);

        if (n > prevn) {
            cout << "n=" << n << endl;
            prevn = n;
            edgeSets.clear();
            excluded[maxdepth] = n + maxdepth;
            genEdgeSets(maxdepth - 1);
        }
        if (!USE_CONCRETE_TREE && (ttt % 100 == 0)) {
            cerr << tt << ",";
        }

        for (size_t v = 0; v < n; ++v) {
            parent[v] = treeDatabase[tt][v + 1];
            children[v].clear();
        }
        for (size_t v = 1; v < n; ++v)
            children[parent[v]].push_back(v);

        type[0] = 0;
        for (size_t v = 1; v < n; ++v)
            type[v] = 1 - type[parent[v]];

        for (size_t v = 0; v < n; ++v)
            gr[v].clear();
        for (size_t v = 1; v < n; ++v) {
            gr[parent[v]].push_back(v);
            gr[v].push_back(static_cast<size_t>(parent[v]));
        }
        for (size_t v = 0; v < n; ++v) {
            deg[v] = gr[v].size();
        }

        res = checkCatterpillarness();
        if (res)
            continue;

        cout << "tree: ";
        for (size_t j = 1; j < n; ++j)
            cout << parent[j] << "->" << j << "; ";
        cout << endl;

        buildLinks();

        bool hadKol0 = false;

        for (size_t setNum = 0; setNum < edgeSets.size(); ++setNum) {
            //----------------- init phase ---------------//
            for (size_t i = 0; i < maxn; ++i)
                needEdge[i] = false; // rename
            m = n;
            if (USE_CONCRETE_TREE) {
                /*for (auto& s: edgeSet) {
                    needEdge[s] = true;
                    m = max(m, s + 1);
                }*/
            } else {
                for (auto& e: edgeSets[setNum]) {
                    needEdge[e] = true;
                    m = max(m, e + 1);
                }
            }

            betaLabels.clear();
            betaSchemes.clear();
            kol = 0; // rename

            cout << "edgeSet={";
            for (size_t e = 0; e < maxn; ++e) {
                if (needEdge[e]) {
                    cout << e << " ";
                }
            }
            cout << "}" << endl;

            for (size_t j = 1; j < n; ++j) {
                if (!isCopy[j]) {
                    for (size_t i = 0; i <= m; ++i) {
                        usedVertex[i] = false;
                        usedEdge[i] = false;
                    }

                    markedVertex[0] = static_cast<size_t>(parent[j]);
                    markedVertex[1] = j;

                    if (type[markedVertex[0]] == 0) {
                        vertexLabels[markedVertex[0]] = 0;
                        vertexLabels[markedVertex[1]] = m - 1;
                    }
                    else {
                        vertexLabels[markedVertex[0]] = m - 1;
                        vertexLabels[markedVertex[1]] = 0;
                    }

                    res = gen(0);
                    if (res) {
                        break;
                    }
                }
            }

            // Conjecture: we should always have at least one k-beta+ labeling, so kol > 0
            if (!kol) {
                if (!hadKol0) {
                    cerr << "OHSHI~" << endl; // rename
                    cout << "tree: ";
                    for (size_t j = 1; j < n; ++j)
                        cout << parent[j] << "->" << j << "; ";
                    cout << endl;
                    countDiam();
                    int maxDeg = 0;
                    for (size_t i = 0; i < n; ++i)
                        maxDeg = max(maxDeg, deg[i]);
                    cout << "diam=" << diam << "; maxDeg=" << maxDeg << endl;
                    cout << "kol=" << kol << ";" << endl;
                }
                hadKol0 = true;
                cout << "edgeSet={";
                for (size_t e = 0; e < maxn; ++e) {
                    if (needEdge[e]) {
                        cout << e << " ";
                    }
                }
                cout << "}" << endl;
            }
        }
        if (hadKol0) {
            cout << endl;
            cout << endl;
            cout << endl;
        }
    }

	endClock = clock();
	double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;
	cerr << "Time: " << elapsed_secs  << "s" << endl;
	cout << "Time: " << elapsed_secs  << "s" << endl;

	return 0;
}

