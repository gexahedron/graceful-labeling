/*
 * File:   beta_plus_seq.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 05 May 2014
 */

// labeling conventions: 
// V \in [0, 1, ..., n-1]
// E \in [1, ..., n-1]
// n < maxn
// We skip caterpillars, but silently assume, 
// that all restrictions work for the most natural alpha-labeling of theirs
// Current set of restrictions:
// beta+
// seq
// We try to restrict with:
// arrowness
// works for n <= 15

// TODO: use conventions for loop variables
// v, u - for vertices; e - for edges

#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <map>
#include "freetree.h"

using namespace std;

// TODO: needs refactoring; put everything into classes and structs and namespaces;
// Or maybe just split everything into blocks commented with their functions-experiments;
const int maxn = 30;
int n;
vector <int> gr[maxn];
int label[maxn];
int seqLabel[maxn];
vector <vector <int> > betaLabels;
vector <vector <int> > betaSeqLabels;
vector <vector <int> > betaSchemes;
vector <int> curBetaLabel;
vector <int> curBetaSeqLabel;
vector <int> curBetaScheme;
int parent[maxn];
vector <int> children[maxn];
pair <int, int> bracketTime[maxn];
pair <int, int> edgeLabels[maxn];
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


int pr[2];
bool res;

int row1[maxn];
int row2[maxn];
int sortedTypes[maxn];
bool halfA, halfB;
int kol;

int curTime;
vector<bool> betaIsArrows;
int kol2;

clock_t beginClock, curClock, prevClock, endClock;

void countDiam() {
	for (int i = 1; i < n; i++)
		if (parent[i] == 0) {
			branch[i] = i;
			diams[i] = 1;
		}
		else
			branch[i] = branch[parent[i]];
	distTo0[0] = 0;
	for (int i = 1; i < n; i++)
		distTo0[i] = distTo0[parent[i]]+1;
	for (int i = 1; i < n; i++)
		diams[branch[i]] = max(diams[branch[i]], distTo0[i]);
	diam = diams[gr[0][0]]+diams[gr[0][1]];
}

bool checkTrees(int v1, int v2) {
	if (deg[v1] != deg[v2])
		return false;
	for (int i = 1; i < (int) gr[v1].size(); i++)
		if (!(checkTrees(gr[v1][i], gr[v2][i])))
			return false;
	return true;
}

void markAsCopy(int v1) {
	isCopy[v1] = true;
	for (int i = 1; i < (int) gr[v1].size(); i++)
		markAsCopy(gr[v1][i]);
}

void buildLinks() {
	for (int i = 0; i < n; i++) {
		link[i] = -1;
		isCopy[i] = false;
	}
	for (int center = 0; center < n; center++) {
		for (int ii = 0; ii < (int) gr[center].size()-1; ii++) {
			int neib1 = gr[center][ii];
			if (neib1 < center)
				continue;
			int neib2 = gr[center][ii+1];
			if (checkTrees(neib1, neib2)) {
				link[neib2] = neib1;
				if (!isCopy[neib2])
					markAsCopy(neib2);
			}
		}
	}

	for (int i = 0; i < n; i++)
		revLinkCount[i] = 0;
	for (int i = n-1; i >= 0; i--)
		if (link[i] != -1)
			revLinkCount[link[i]] = revLinkCount[i]+1;
}

// TODO
bool checkAlpha() {
  return false;
}

bool checkSequentiallity() {
    // TODO: split checkSequentiallity into checkAlpha and checkBetaSeq
    // And always return true if good tree, this true would be analyzed further in gen()

	/*int maxT0, maxT1;
	int minT0, minT1;
	maxT0 = maxT1 = 0;
	minT0 = minT1 = 2*n;
	for (int i = 0; i < n; i++)
		if (type[i] == 0) {
			maxT0 = max(maxT0, label[i]);
			minT0 = min(minT0, label[i]);
		}
		else {
			maxT1 = max(maxT1, label[i]);
			minT1 = min(minT1, label[i]);
		}

	if (min(maxT0, maxT1) < max(minT0, minT1)) {
		hasAlpha = true;
		kol++;
		return true;
	}*/


	for (int i = 0; i < n; i++)
		if (type[i] == 0) {
			row1[label[i]] = 1;
			row2[n-1 - label[i]] = 0;
		}
		else {
			row1[label[i]] = 0;
			row2[n-1 - label[i]] = 1;
		}

	for (int ii = n-1; ii > 0; ii--) {
		int countOverlaps = 0;
		int whereOverlap;
		for (int j = 0; j < n-1; j++)
			if (row1[j] && row2[j]) {
				countOverlaps++;
				whereOverlap = j;
				if (countOverlaps > 1)
					break;
			}
		if (countOverlaps == 1) {
			//it is harmonious, now let's check sequentiallity
			
			//use whereOverlap
			//cout << endl;
			//cout << ii << " -- " << whereOverlap << endl;
			int tp0, tp1;
			for (int i = 0; i < n; i++)
				if (type[i] == 0) {
					seqLabel[i] = (label[i] + n-1-whereOverlap) % (n-1);
					if (seqLabel[i] == 0)
						tp0 = i;
				}
				else {
					seqLabel[i] = (n-1-label[i] + n-1-whereOverlap + ii) % (n-1);
					if (seqLabel[i] == 0)
						tp1 = i;
				}

			
			bool isSeq0, isSeq1;
			isSeq0 = false;
			isSeq1 = false;
			int minv, maxv;
			seqLabel[tp0] = n-1;
			seqLabel[tp1] = 0;
			minv = seqLabel[1]+seqLabel[parent[1]];
			maxv = minv;
			for (int i = 2; i < n; i++) {
				int e = seqLabel[i] + seqLabel[parent[i]];
				minv = min(minv, e);
				maxv = max(maxv, e);
			}
			isSeq0 = (maxv-minv+1 == n-1);

			
			seqLabel[tp0] = 0;
			seqLabel[tp1] = n-1;
			minv = seqLabel[1]+seqLabel[parent[1]];
			maxv = minv;
			for (int i = 2; i < n; i++) {
				int e = seqLabel[i] + seqLabel[parent[i]];
				minv = min(minv, e);
				maxv = max(maxv, e);
			}
			isSeq1 = (maxv-minv+1 == n-1);
			
			if (isSeq0 || isSeq1) {
				if (isSeq0) {
					seqLabel[tp0] = n-1;
					seqLabel[tp1] = 0;
				}
				bool half0 = true;
				bool half1 = true;
				bool half0Inv = true;
				bool half1Inv = true;
				for (int i = 0; i < n; i++)
					if (type[i] == 0) {
						half0 = half0 && (label[i] == seqLabel[i]);
						half0Inv = half0Inv && (label[i] == (n-1-seqLabel[i]));
					}
					else {
						half1 = half1 && (label[i] == seqLabel[i]);
						half1Inv = half1Inv && (label[i] == (n-1-seqLabel[i]));
					}

				if (!half0 && !half0Inv && !half1 && !half1Inv && isSeq1) {
					seqLabel[tp0] = 0;
					seqLabel[tp1] = n-1;

					half0 = true;
					half1 = true;
					half0Inv = true;
					half1Inv = true;
					for (int i = 0; i < n; i++)
						if (type[i] == 0) {
							half0 = half0 && (label[i] == seqLabel[i]);
							half0Inv = half0Inv && (label[i] == (n-1-seqLabel[i]));
						}
						else {
							half1 = half1 && (label[i] == seqLabel[i]);
							half1Inv = half1Inv && (label[i] == (n-1-seqLabel[i]));
						}
				}

				if (!half0 && !half0Inv && !half1 && !half1Inv)
					return false;

				if (half0 || half0Inv)
					halfA = true;
				if (half1 || half1Inv)
					halfB = true;
				kol++;

				if (half0Inv || half1Inv) {
					for (int i = 0; i < n; i++)
						seqLabel[i] = n-1-seqLabel[i];
				}

				curBetaLabel.clear();
				curBetaSeqLabel.clear();
				curBetaScheme.clear();

				for (int i = 0; i < n; i++) {
					curBetaLabel.push_back(label[i]);
					curBetaSeqLabel.push_back(seqLabel[i]);
					curBetaScheme.push_back(0);
				}
				for (int i = 0; i < n; i++) {
					curBetaScheme[label[i]] = type[i];
				}
				return true;
			}
		}
		int tt = row2[0];
		for (int j = 0; j < n-2; j++)
			row2[j] = row2[j+1];
		row2[n-2] = tt;
	}
	return false;
}

bool checkArrowness() {
    for (int v = 1; v < n; ++v) {
        int e = abs(curBetaLabel[v] - curBetaLabel[parent[v]]);
        edgeLabels[e].first = parent[v];
        edgeLabels[e].second = v;
    }
    bool hasOutgoingArrow[maxn];
    for (int v = 0; v < n; ++v)
        hasOutgoingArrow[v] = false;
    for (int e1 = 1; e1 < n - e1; ++e1) {
        int e2 = n - e1;
        int v11 = edgeLabels[e1].first;
        int v12 = edgeLabels[e1].second;
        int v21 = edgeLabels[e2].first;
        int v22 = edgeLabels[e2].second;
        int ar1, ar2;
        int time1 = bracketTime[v11].first;
        int time2 = bracketTime[v21].first;
        if (time1 < time2) {
            ar1 = v22;
            if (bracketTime[v12].first <= time2 && bracketTime[v21].second <= bracketTime[v12].second)
                ar2 = v11;
            else
                ar2 = v12;
        } else {
            ar1 = v12;
            if (bracketTime[v22].first <= time1 && bracketTime[v11].second <= bracketTime[v22].second)
                ar2 = v21;
            else
                ar2 = v22;
        }
        //cout << ar1 << " " << ar2 << endl;
        if (hasOutgoingArrow[ar1])
            return false;
        hasOutgoingArrow[ar1] = true;
        if (hasOutgoingArrow[ar2])
            return false;
        hasOutgoingArrow[ar2] = true;
    }
    if (!(n % 2)) {
        int e = n / 2;
        int v1 = edgeLabels[e].first;
        int v2 = edgeLabels[e].second;
        //cout << v1 << " " << v2 << endl;
        if (hasOutgoingArrow[v1] && hasOutgoingArrow[v2])
            return false;
    }
    //cout << endl;
    return true;
}

bool gen(int v) {
	if (v == n) {
        // TODO: check isAlpha
        bool isSeq = checkSequentiallity();
        if (isSeq) {
            betaLabels.push_back(curBetaLabel);
			betaSeqLabels.push_back(curBetaSeqLabel);
			betaSchemes.push_back(curBetaScheme);
            bool isArr = checkArrowness();
            if (isArr)
                ++kol2;
            betaIsArrows.push_back(isArr);
        }
        return false;
        //return (isSeq || isAlpha);
	}
	else {
		int edgeVal;
		int lowVal = 1;
		int hiVal = n-2;
		int startVal = lowVal;
		int endVal = hiVal;
		int dir = 1;

		if ((v != pr[0]) && (v != pr[1])) {
			if (type[v] == 0) {
				if (parent[v] != -1)
					hiVal = label[parent[v]]-1;
				startVal = lowVal;
				endVal = hiVal;
				dir = 1;
				if (startVal > endVal)
					return false;
			}
			else {
				if (parent[v] != -1)
					lowVal = label[parent[v]]+1;
				startVal = hiVal;
				endVal = lowVal;
				dir = -1;
				if (startVal < endVal)
					return false;
			}
		}


		for (int i = startVal; i != endVal+dir; i+=dir) {
			if ((usedVertex[i]) && (v != pr[0]) && (v != pr[1]))
				continue;
			else {
				if ((v != pr[0]) && (v != pr[1]))
					label[v] = i;

				if (parent[v] != -1) {
					edgeVal = abs(label[v] - label[parent[v]]);
					if (usedEdge[edgeVal]) {
						if ((v == pr[0]) || (v == pr[1])) {
							break;
					  } else {
							continue;
            }
          }

					if (link[v] != -1) {
						if (edgeVal > abs(label[link[v]]-label[parent[link[v]]])) {
							continue;
            }
          }
					if (edgeVal < revLinkCount[v]) {
						continue;
          }
					usedEdge[edgeVal] = true;
				}
				usedVertex[label[v]] = true;

				if (gen(v+1))
					return true;

				if (parent[v] != -1)
					usedEdge[edgeVal] = false;

				usedVertex[label[v]] = false;
				if ((v == pr[0]) || (v == pr[1]))
					break;
			}
		}
	}
	return false;
}


bool checkCatterpillarness() {
	int curV, nextV, prevV;
	for (int v = 0; v < n; v++) 
		if (deg[v] > 1) {
			int kolDegNot1 = 0;
			for (int i = 0; i < (int) gr[v].size(); i++) {
				if (deg[gr[v][i]] > 1) {
					kolDegNot1++;
					nextV = gr[v][i];
				}
			}
			if (kolDegNot1 == 1) {
				curV = v;
				break;
			}
			if (kolDegNot1 == 0)
				return true;
		}
	while (curV != nextV) {
		prevV = curV;
		curV = nextV;
		int kolDegNot1 = 0;
		for (int i = 0; i < (int) gr[curV].size(); i++) {
			int u = gr[curV][i];
			if ((deg[u] > 1) && (u != prevV)) {
				nextV = u;
				kolDegNot1++;
			}
		}
		if (kolDegNot1 > 1)
			return false;
	}
	return true;
}

void fillBracketTimes(int v) {
    bracketTime[v].first = curTime;
    ++curTime;
    for (int i = 0; i < children[v].size(); ++i)
        fillBracketTimes(children[v][i]);
    bracketTime[v].second = curTime;
    ++curTime;
}

int main(int argc, char* argv[])
{
	//bool tempCode = freopen("maps7-14_verbose.txt", "w", stdout);
	bool tempCode = freopen("test.txt", "w", stdout);
	beginClock = clock();
	prevClock = beginClock;
	srand (time(NULL));

	int razmBegin = 7;
	int razmEnd = 9;
	for (int i = razmBegin; i <= razmEnd; i++)
		freeTree(i);
	cerr << "#trees=" << totalAmount << endl;
	
    /*
    //////////////////// testing concrete trees here /////////////////////
    totalAmount = 1;
    //const int N = 16;
    //int curTree[N+1] = {N, -1, 0, 1, 2, 3, 4, 0, 6, 7, 8, 9, 0, 11, 12, 13, 14};
    const int N = 25;
    int curTree[N+1] = {N, -1, 0, 1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0, 21, 0, 23};
    for (size_t i = 0; i < N+1; ++i)
        treeDatabase[0][i] = curTree[i];
    //////////////////////////////////////////////////////////////////////
    */

	int prevn = 0;
	for (int ttt = 0; ttt < totalAmount; ttt++) {
		int tt = ttt;
		n = treeDatabase[tt][0];
		
		if (n > prevn) {
			cout << "n=" << n << endl;
			prevn = n;
		}
		
		if (ttt % 100 == 0)
			cerr << tt << ",";

		for (int v = 0; v < n; ++v)
			parent[v] = treeDatabase[tt][v+1];

        for (int v = 0; v < n; ++v)
            children[v].clear();
        for (int v = 1; v < n; ++v)
            children[parent[v]].push_back(v);
        curTime = 0;
        fillBracketTimes(0);
		type[0] = 0;
		for (int i = 1; i < n; i++)
			type[i] = 1 - type[parent[i]];

		for (int i = 0; i < n; i++)
			gr[i].clear();

		for (int i = 1; i < n; i++) {
			gr[parent[i]].push_back(i);
			gr[i].push_back(parent[i]);
		}

		for (int i = 0; i < n; i++)
			deg[i] = 0;

		for (int i = 1; i < n; i++) {
			deg[parent[i]]++;
			deg[i]++;
		}

		betaLabels.clear();
		betaSeqLabels.clear();
		betaSchemes.clear();

        betaIsArrows.clear();

		res = checkCatterpillarness();
		if (!res) {
			buildLinks();
			halfA = false;
			halfB = false;
			kol = 0;
            kol2 = 0;
			hasAlpha = false;
			
			for (int j = 1; j < n; j++) {
				if (!isCopy[j]) {
					for (int i = 0; i <= n; i++) {
						usedVertex[i] = false;
						usedEdge[i] = false;
					}
					usedEdge[0] = true;

					pr[0] = parent[j];
					pr[1] = j;

					if (type[pr[0]] == 0) {
						label[pr[0]] = 0;
						label[pr[1]] = n-1;
					} 
					else {
						label[pr[0]] = n-1;
						label[pr[1]] = 0;
					}
			
					res = gen(0);
					if (hasAlpha)
						break;
				}
            }

				if (!kol) { // We always have at least 1 beta-seq => kol > 0 always
						cout << "OHSHI~" << endl;
						cerr << "OHSHI~" << endl;
				}
				if (kol) {
					// print out all labelings
                    for (int i = 0; i < betaLabels.size(); i++) {
						cout << "beta labels: ";
						for (int j = 0; j < n; j++)
							cout << betaLabels[i][j] << " ";
						cout << endl;
						
                        cout << "seq labels: ";
						for (int j = 0; j < n; j++)
							cout << betaSeqLabels[i][j] << " ";
						cout << endl;
						
                        cout << "scheme: ";
						for (int j = 0; j < n; j++)
							cout << betaSchemes[i][j] << " ";
						cout << endl;
						
                        cout << "arrowness: " << betaIsArrows[i] << endl;

                        cout << "|";
						cout << endl;
					}

					cout << "tree: ";
					for (int j = 1; j < n; j++)
						cout << parent[j] << "->" << j << "; ";
					cout << endl;
					countDiam();
					int maxDeg = 0;
					for (int i = 0; i < n; i++)
						maxDeg = max(maxDeg, deg[i]);
					cout << "diam=" << diam << "; maxDeg=" << maxDeg << endl;
                    cout << "kol=" << kol << "=>" << kol2 << ";" << endl;

					int sizeL, sizeU;
					sizeL = 0;
					sizeU = 0;
					for (int i = 0; i < n; i++)
						if (betaSchemes[0][i] == 0)
							sizeL++;
						else
							sizeU++;
					int alphaSum = n*(n-1)/2 - sizeL*(sizeL-1);
                    
                    int maxAlphaDiff = 0;
                    int restrictedMaxAlphaDiff = 0;

                    int minAlphaDiff = maxn*maxn*2;
                    int restrictedMinAlphaDiff = minAlphaDiff;

                    map<int,int> diffs;
                    map<int,int> restrictedDiffs;

					for (int i = 0; i < betaSchemes.size(); i++) {
						cout << "scheme: ";
						int curSum = 0;
						for (int j = 0; j < n; j++) {
							cout << betaSchemes[i][j] << " ";
							if (betaSchemes[i][j] == 0)
								curSum -= j;
							else
								curSum += j;
						}
                        int alphaDiff = (alphaSum - curSum) / 2;
                        
                        maxAlphaDiff = max(maxAlphaDiff, alphaDiff);
                        minAlphaDiff = min(minAlphaDiff, alphaDiff);
                        ++diffs[alphaDiff];

                        if (betaIsArrows[i]) {
                            restrictedMaxAlphaDiff = max(restrictedMaxAlphaDiff, alphaDiff);
                            restrictedMinAlphaDiff = min(restrictedMinAlphaDiff, alphaDiff);
                            ++restrictedDiffs[alphaDiff];
                        }

						cout << "; ";
						cout << "sizes: " << sizeL << " " << sizeU << "; ";
						cout << "alphaDiff=" << alphaDiff;
						cout << "; ";
                        cout << "arrowness=" << betaIsArrows[i];
                        cout << "; ";
						cout << endl;
					}

                    cout << "perc=" << 1.0*kol2/kol << endl;
                    if (!kol2)
                        cout << "Achievement! FAIL" << endl;
                    if (!maxAlphaDiff)
                        cout << "Achievement! ALPHA ONLY" << endl;
                    if (kol2 && maxAlphaDiff && !restrictedMaxAlphaDiff) {
                        cout << "Achievement! NEW ALPHA" << endl;
                        if (kol2 % 2)
                            cout << "Achievement! BUGS" << endl;
                    }
                    if (kol == 1 || (kol == 2 && !maxAlphaDiff))
                        cout << "Achievement! UNIQUE ONLY" << endl;
                    if (kol2 < kol && (kol2 == 1 || (kol2 == 2 && !restrictedMaxAlphaDiff)))
                        cout << "Achievement! NEW UNIQUE" << endl;

                    if (minAlphaDiff == maxAlphaDiff && minAlphaDiff > 0)
                        cout << "Achievement! TYPED ONLY: " << maxAlphaDiff << endl;
                    if (kol2 && minAlphaDiff < restrictedMinAlphaDiff)
                        cout << "Achievement! KILLED MIN: " << minAlphaDiff << "~>" << restrictedMinAlphaDiff << endl;
                    if (kol2 && maxAlphaDiff > restrictedMaxAlphaDiff)
                        cout << "Achievement! KILLED MAX: " << maxAlphaDiff << "~>" << restrictedMaxAlphaDiff << endl;
                    if (kol2 && !minAlphaDiff && restrictedMinAlphaDiff)
                        cout << "Achievement! KILLED ALPHA" << endl;
                    if (kol2 && minAlphaDiff != maxAlphaDiff && restrictedMinAlphaDiff == restrictedMaxAlphaDiff)
                        cout << "Achievement! NEW TYPED: " << "("<< minAlphaDiff << ", " << maxAlphaDiff << ")~>" << restrictedMinAlphaDiff << endl;
                    if (diffs.size() != restrictedDiffs.size()) {
                        cout << "Achievement! KILLED SOME DIFFS: " << "[";
                        for (map<int,int>::const_iterator iter = diffs.begin(); iter != diffs.end(); ++iter) {
                            if (iter != diffs.begin())
                                cout << ",";
                            cout << "{" << iter->first << "; " << iter->second << "}";
                        }

                        cout << "] ~> [";
                        for (map<int,int>::const_iterator iter = restrictedDiffs.begin(); iter != restrictedDiffs.end(); ++iter) {
                            if (iter != restrictedDiffs.begin())
                                cout << ",";
                            cout << "{" << iter->first << "; " << iter->second << "}";
                        }
                        cout << "]" << endl;
                    }
					cout << endl;
					cout << endl;
					cout << endl;
				}
		}
	}

	endClock = clock();
	double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;

	cerr << "Time: " << elapsed_secs  << "s" << endl;
	cout << "Time: " << elapsed_secs  << "s" << endl;
	cerr << "the end" << endl;
	//system ("pause");
	return 0;
}

