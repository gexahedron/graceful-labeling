/*
 * File:   edge_rotatability.cpp
 * Author: Nikolay Ulyanov
 *
 * Created somewhere in Autumn 2013; and again for seq-labeling on 09 July 2015
 */


#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include "freetree.h"

using namespace std;

const int maxn = 30;
int n;
vector <int> gr[maxn];
int label[maxn];
int parent[maxn];
bool usedVertex[maxn];
bool usedEdge[maxn];
short treeDatabase[30000000][30];
vector <int> newTree;
int totalAmount;

int pr[2];
bool res[maxn];

bool gen(int v) {
	//cout << v << endl;
	if (v == n)
		return true;
	else {
		int edgeVal;
		for (int i = n-1; i > 1; i--) {
			if ((usedVertex[i]) && (v != pr[0]) && (v != pr[1]))
				continue;
			else {
				if ((v != pr[0]) && (v != pr[1]))
					label[v] = i;

				if (parent[v] != -1) {
					edgeVal = abs(label[v] - label[parent[v]]);
					if (usedEdge[edgeVal])
						if ((v == pr[0]) || (v == pr[1]))
							break;
						else
							continue;
					usedEdge[edgeVal] = true;
				}
				usedVertex[label[v]] = true;

				/*for (int j = 0; j < n; j++)
					cout << label[j] << " ";
				cout << endl;*/

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

int main(int argc, char* argv[])
{
	clock_t begin = clock();
	srand (time(NULL));
	totalAmount = 0;
	
	int sdvig = 0;

	cout << "graph { " << endl;
	cout << "    concentrate=true;" << endl;
	cout << "    overlap=true;" << endl;
	cout << "    rankdir=LR;" << endl;

    int razm = 8;
	cin >> razm;
	for (int i = razm; i <= razm; i++)
		freeTree(i);
	cerr << "#trees=" << totalAmount << endl;
	
	int prevn = 0;
	for (int ttt = 0; ttt < totalAmount; ttt++) {
		
		int tt = ttt;
		n = treeDatabase[tt][0];
		
		if (n > prevn) {
			//cout << "n=" << n << endl;
			prevn = n;
		}
		
		if (ttt % 100 == 0)
			cerr << tt << ",";

		for (int i = 0; i < n; i++)
			gr[i].clear();
		
		for (int i = 0; i < n; i++)
			parent[i] = treeDatabase[tt][i+1];

		for (int i = 1; i < n; i++) {
			gr[parent[i]].push_back(i);
		}
		
		usedEdge[0] = true;

		int kol = 0;
		for (int j = 1; j < n; j++) {
			for (int i = 1; i <= n; i++) {
				usedVertex[i] = false;
				usedEdge[i] = false;
			}

			pr[0] = parent[j];
			pr[1] = j;
			label[pr[0]] = 1;
			label[pr[1]] = n;
			
			res[j] = gen(0);
			if (!res[j])
				kol++;
			/*if (!gen(0)) {
				cout << tt << ": ";
				cout << "(" << pr[0] << ", " << pr[1] << ") ";
				for (int i = 2; i <= n; i++)
					cout << i-1 << "->" << treeDatabase[tt][i] << "; ";
				cout << endl;
			}*/
			//cout << endl;
		}

		if (kol > 0) {
			for (int j = 1; j < n; j++) {
				cout << "    " << parent[j] + sdvig << " -- " << j + sdvig;
				if (!res[j])
					cout << "[style=dashed,color=red,penwidth=3.0]";
				cout << ";" << endl;
			}
			cout << endl;
			sdvig += n;
		}
	}
	cout << "}" << endl;
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cerr << "Time: " << elapsed_secs  << "s" << endl;
	cerr << "the end" << endl;
	return 0;
}

