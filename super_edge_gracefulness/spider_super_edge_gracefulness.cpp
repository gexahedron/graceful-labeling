/*
 * File:   spider_super_edge_gracefulness.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 14 October 2016
 */

#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <set>
#include <cassert>
#include <algorithm>

using namespace std;

const int maxn = 30;

short treeDatabase[10000][maxn];
short treeLengths[10000][maxn];
short curTreeLengths[maxn];
static int totalAmount = 0;
int lengths[maxn];

int n;
vector<int> gr[maxn];
int label[maxn];
int order[maxn];
vector<vector<int>> edgeGracefulLabels;
vector<int> curEdgeGracefulLabel;

int parent[maxn];
vector<int> children[maxn];
pair<int, int> edgeLabels[maxn];
//int type[maxn];
int deg[maxn];
int leaf_count[maxn];
int link[maxn];
bool isCopy[maxn];
int revLinkCount[maxn];
bool usedVal[maxn];
int edges[maxn];
int cur_sum[maxn];

int diam;
int branch[maxn];
int distTo0[maxn];
int diams[maxn];

bool used_e_at_v[maxn][maxn];
int count_used[maxn];
bool res;
int pairs[maxn][2];
vector<string> profiles;
set<string> all_profiles;
set<string> solved_profiles;
map<string, int> prof_count;
vector<string> schemes;
map<string, int> scheme_count;
int sol_count;
string profile;
set<string> sols;
bool all_lengths_different;

clock_t beginClock, curClock, prevClock, endClock;


void countDiam() {
  for (int i = 1; i < n; ++i) {
    if (parent[i] == 0) {
      branch[i] = i;
      diams[i] = 1;
    } else {
      branch[i] = branch[parent[i]];
    }
  }
  distTo0[0] = 0;
  for (int i = 1; i < n; ++i) {
    distTo0[i] = distTo0[parent[i]] + 1;
  }
  for (int i = 1; i < n; ++i) {
    diams[branch[i]] = max(diams[branch[i]], distTo0[i]);
  }
  diam = diams[gr[0][0]] + diams[gr[0][1]];
}

bool checkTrees(int v1, int v2) {
  if (deg[v1] != deg[v2]) {
    return false;
  }
  for (int i = 1; i < (int) gr[v1].size(); ++i) {
    if (!(checkTrees(gr[v1][i], gr[v2][i]))) {
      return false;
    }
  }
  return true;
}

void markAsCopy(int v1) {
  isCopy[v1] = true;
  for (int i = 1; i < (int) gr[v1].size(); ++i) {
    markAsCopy(gr[v1][i]);
  }
}

void buildLinks() {
  for (int i = 0; i < n; ++i) {
    link[i] = -1;
    isCopy[i] = false;
  }
  for (int center = 0; center < n; ++center) {
    for (int ii = 0; ii < (int) gr[center].size() - 1; ++ii) {
      int neib1 = gr[center][ii];
      if (neib1 < center) {
        continue;
      }
      int neib2 = gr[center][ii + 1];
      if (checkTrees(neib1, neib2)) {
        link[neib2] = neib1;
        if (!isCopy[neib2]) {
          markAsCopy(neib2);
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    revLinkCount[i] = 0;
  }
  for (int i = n - 1; i >= 0; --i) {
    if (link[i] != -1) {
      revLinkCount[link[i]] = revLinkCount[i] + 1;
    }
  }
}

int segabs(int v) {
  if (v >= 0) {
    return v;
  }
  return v + n;
}

bool gen(int e) {
  if (e == n) {
    vector<string> cur_sol1_parts, cur_sol2_parts;
    string part1, part2;
    for (int i = 1; i < n; ++i) {
      part1 += to_string(label[i]) + ",";
      part2 += to_string(-label[i]) + ",";
      if (deg[i] == 1) {
        part1 += ";";
        part2 += ";";
        cur_sol1_parts.push_back(part1);
        cur_sol2_parts.push_back(part2);
        part1 = "";
        part2 = "";
      }
    }
    sort(cur_sol1_parts.begin(), cur_sol1_parts.end());
    sort(cur_sol2_parts.begin(), cur_sol2_parts.end());
    string cur_sol1, cur_sol2;
    for (const auto& p : cur_sol1_parts) {
      cur_sol1 += p;
    }
    for (const auto& p : cur_sol2_parts) {
      cur_sol2 += p;
    }

    if (sols.find(cur_sol1) != sols.end()) {
      return false;
    }
    if (sols.find(cur_sol2) != sols.end()) {
      return false;
    }

    /*//about +- changes in scheme
    //one counterexample: 1,4;1,4
    //but it doesn't work with other heuristics
    int vtypes[maxn];
    for (int i = 1; i < n; ++i) {
      vtypes[i] = 0;
    }
    vtypes[n - 1] = 1;
    while (true) {
      int idx = n - 2;
      bool found = false;
      while (idx > 0) {
        if (deg[idx] != 1 && vtypes[idx] == 0 && vtypes[idx + 1] != 0) {
          found = true;
          vtypes[idx] = -vtypes[idx + 1];
          for (int i = 1; i < n; ++i) {
            if (label[i] == -label[idx]) {
              if (vtypes[i] != 0 && vtypes[i] != -vtypes[idx]) {
                return false;
              }
              vtypes[i] = -vtypes[idx];
              break;
            }
          }
        } else if (deg[idx] != 1 && vtypes[idx] != 0 && vtypes[idx + 1] == 0) {
          found = true;
          vtypes[idx + 1] = -vtypes[idx];
          for (int i = 1; i < n; ++i) {
            if (label[i] == -label[idx + 1]) {
              if (vtypes[i] != 0 && vtypes[i] != -vtypes[idx + 1]) {
                return false;
              }
              vtypes[i] = -vtypes[idx + 1];
              break;
            }
          }
        }
        --idx;
      }
      if (!found) {
        break;
      }
    }
    for (int i = 1; i < n - 1; ++i)
      if (vtypes[i] != 0 && vtypes[i + 1] != 0 && deg[i] != 1 && vtypes[i] == vtypes[i + 1]) {
        return false;
      }*/



    //arithmetic labelings
    //family of counterexamples:
    //1,2,1,2,1,2,1,2;
    //1,2,1,2,1,2,1,2,1,2,1,2;
    //...
    /*map<string, int> triples;
    for (int i = 1; i < n; ++i) {
      if (deg[i] != 1) {
        vector<string> terms;
        terms.push_back(to_string(abs(label[i])));
        terms.push_back(to_string(abs(label[i + 1])));
        terms.push_back(to_string(abs(label[i] + label[i + 1])));
        sort(terms.begin(), terms.end());
        string triple = terms[0] + "," + terms[1] + "," + terms[2];
        ++triples[triple];
        if (triples[triple] > 2) {
          return false;
        }
      }
    }
    if (2 * triples.size() != n - 1 - deg[0]) {
      return false;
    }
*/
    //unique edges on path
    set<int> edges_on_path;
    for (int i = 1; i < n; ++i) {
      if (deg[i] == 1) {
        if (edges_on_path.find(-label[i]) != edges_on_path.end()) {
          return false;
        }

        edges_on_path.clear();
      } else {
        if (edges_on_path.find(-label[i]) != edges_on_path.end()) {
          return false;
        }
        edges_on_path.insert(label[i]);
      }
    }

    /*//edges always follow same order; transitivity
    bool hier[maxn][maxn];
    for (int i = 1; i < n; ++i)
      for (int j = 1; j < n; ++j)
        hier[i][j] = false;
    vector<int> path;
    for (int i = 1; i < n; ++i) {
      path.push_back(abs(label[i]));
      if (deg[i] == 1) {
        for (int j = 0; j < path.size(); ++j) {
          for (int k = j + 1; k < path.size(); ++k) {
            if (hier[path[k]][path[j]]) {
              return false;
            }
          }
        }

        for (int j = 0; j < path.size(); ++j) {
          for (int k = j + 1; k < path.size(); ++k) {
            hier[path[j]][path[k]] = true;
          }
        }
        path.clear();
      }
    }*/


    // for (int i = 1; i < n; ++i) {
    //   if (abs(label[i]) == n / 2 && deg[i] != 1) {
    //     return false;
    //   }
    // }


    if (all_lengths_different) {
      string all_letters = "abcdefghijklmnopqrstuvwxyz";
      map<int, char> letters;
      string cur_scheme = profile + " => ";
      int j = -1;
      for (int ii = 1; ii < n; ++ii) {
        int i = label[ii];
        if (letters.find(-i) != letters.end()) {
          cur_scheme += "-";
          cur_scheme += letters[-i];
        } else {
          ++j;
          letters[i] = all_letters[j];
          cur_scheme += letters[i];
        }

        if (deg[ii] == 1) {
          cur_scheme += "; ";
        } else {
          cur_scheme += " ";
        }
      }
      schemes.push_back(cur_scheme);
      ++scheme_count[cur_scheme];
    }

    sols.insert(cur_sol1);
    sols.insert(cur_sol2);
    solved_profiles.insert(profile);
    ++sol_count;
    ++prof_count[profile];
    curEdgeGracefulLabel.clear();
    for (int i = 1; i < n; ++i) {
      curEdgeGracefulLabel.push_back(label[i]);
    }
    edgeGracefulLabels.push_back(curEdgeGracefulLabel);
    profiles.push_back(profile);
    return false;
  } else {
    int upper = n - e;
    for (int ii = 0; ii < upper; ++ii) {
      int i = edges[ii];
      int v = order[e - 1];
      int next_v = order[e];
      label[v] = i;

      if (link[v] != -1) {
        if (segabs(i) > segabs(label[link[v]])) {
          continue;
        }
      }
      if (segabs(i) < revLinkCount[v]) {
        continue;
      }

      if (deg[v] == 1) {
        if (abs(cur_sum[v] + i) > n / 2) {
          continue;
        }
        if (usedVal[segabs(cur_sum[v] + i)]) {
          continue;
        }
      }

      if (deg[parent[v]] > 2) {
        int add = 0;
        if (!used_e_at_v[parent[v]][abs(i)]) {
          add = 1;
        }
        if (count_used[parent[v]] + add > deg[parent[v]] / 2) {
          continue;
        }

        if (!used_e_at_v[parent[v]][abs(i)]) {
          pairs[abs(i)][0] = curTreeLengths[v];
        } else {
          pairs[abs(i)][1] = curTreeLengths[v];
          if (pairs[abs(i)][0] == 1 && pairs[abs(i)][1] == 1) {
            continue;
          }
          /*if (i > 0 && pairs[abs(i)][0] > pairs[abs(i)][1]) {
            continue;
          }
          if (i < 0 && pairs[abs(i)][0] < pairs[abs(i)][1]) {
            continue;
          }*/
        }

      }

      if (e == n - 1 || parent[next_v] != parent[v]) {
        if (abs(cur_sum[parent[v]] + i) > n / 2) {
          continue;
        }
        if (usedVal[segabs(cur_sum[parent[v]] + i)]) {
          continue;
        }

        if (deg[parent[v]] > 2) {
          if (segabs(cur_sum[parent[v]] + i) != 0) {
            continue;
          }
          int maxsumlen = 0;
          bool has_half = false;
          for (int ee = 0; ee <= n / 2; ++ee) {
            if (used_e_at_v[parent[v]][ee]) {
              if (pairs[ee][0] + pairs[ee][1] >= maxsumlen) {
                if (pairs[ee][0] + pairs[ee][1] > maxsumlen) {
                  has_half = false;
                }
                maxsumlen = pairs[ee][0] + pairs[ee][1];
                if (abs(pairs[ee][0] - pairs[ee][1]) <= 1) {
                  has_half = true;
                }
              }
            }
          }
          if (!has_half) {
            continue;
          }

          vector<string> profile_parts;
          for (int ee = 0; ee <= n / 2; ++ee) {
            if (used_e_at_v[parent[v]][ee]) {
              profile_parts.push_back(to_string(min(pairs[ee][0], pairs[ee][1])) + "," +
                  to_string(max(pairs[ee][0], pairs[ee][1])));
            }
          }
          sort(profile_parts.begin(), profile_parts.end());
          profile = "";
          for (int j = 0; j < profile_parts.size(); ++j)
            profile += profile_parts[j] + ";";
          // if (solved_profiles.find(profile) != solved_profiles.end()) {
          //   continue;
          // }
          all_profiles.insert(profile);

          // if (!used_e_at_v[parent[v]][1]) {
          //   continue;
          // }

          vector<int> root_edges;
          vector<int> sums;
          for (int ee = 0; ee <= n / 2; ++ee) {
            if (used_e_at_v[parent[v]][ee]) {
              sums.push_back(pairs[ee][0] + pairs[ee][1]);
              root_edges.push_back(ee);
            }
          }
          sort(root_edges.begin(), root_edges.end(), std::greater<int>());
          sort(sums.begin(), sums.end());
          int cur_sum = 0;
          bool specific_edges = true;
          for (int i = 0; i < sums.size(); ++i) {
            cur_sum += sums[i];
            if (root_edges[i] != (n - cur_sum + 1) / 2) {
              specific_edges = false;
              break;
            }
          }
          // if (!specific_edges) {
          //   continue;
          // }
        } else if (segabs(cur_sum[parent[v]] + i) == 0) {
          continue;
        }
      }

      cur_sum[v] += i;
      cur_sum[parent[v]] += i;
      if (deg[v] == 1) {
        usedVal[segabs(cur_sum[v])] = true;
      }

      bool added_par = false;
      if (!used_e_at_v[parent[v]][abs(i)]) {
        added_par = true;
        count_used[parent[v]] += 1;
        used_e_at_v[parent[v]][abs(i)] = true;
      }

      bool added = false;
      if (!used_e_at_v[v][abs(i)]) {
        added = true;
        count_used[v] += 1;
        used_e_at_v[v][abs(i)] = true;
      }

      if (e == n - 1 || parent[next_v] != parent[v]) {
        usedVal[segabs(cur_sum[parent[v]])] = true;
      }
      swap(edges[ii], edges[n - e - 1]);

      if (gen(e + 1)) {
        return true;
      }

      swap(edges[ii], edges[n - e - 1]);
      if (deg[v] == 1) {
        usedVal[segabs(cur_sum[v])] = false;
      }
      if (e == n - 1 || parent[next_v] != parent[v]) {
        usedVal[segabs(cur_sum[parent[v]])] = false;
      }

      if (added_par) {
        count_used[parent[v]] -= 1;
        used_e_at_v[parent[v]][abs(i)] = false;
      }
      if (added) {
        count_used[v] -= 1;
        used_e_at_v[v][abs(i)] = false;
      }

      cur_sum[v] -= i;
      cur_sum[parent[v]] -= i;

      /*if (solved_profiles.find(profile) != solved_profiles.end()) {
        if (e > deg[0]) {
          break;
        } else {
          continue;
        }
      }*/

    }
  }
  return false;
}

void genTree(int iter, int num, int min_len) {
  if (num == 0) {
    if (iter < 3 || iter % 2 != 0) {
      return;
    }
    /*if (lengths[iter - 2] >= n / 2 - 1) {
      return;
    }*/
    /*if (lengths[iter - 1] > n / 2 - 1) {
      return;
    }*/
    treeDatabase[totalAmount][0] = n;
    treeDatabase[totalAmount][1] = -1;
    int cur_count = 1;
    int one_count = 0;
    for (int i = 0; i < iter; ++i) {
      if (lengths[i] == 1) {
        ++one_count;
      }
      treeDatabase[totalAmount][cur_count + 1] = 0;
      treeLengths[totalAmount][cur_count] = lengths[i];
      ++cur_count;
      for (int j = 0; j < lengths[i] - 1; ++j) {
        treeDatabase[totalAmount][cur_count + 1] = cur_count - 1;
        ++cur_count;
      }
    }
    if (one_count > iter / 2) {
      return;
    }
    ++totalAmount;
    return;
  }
  for (int i = min_len; i <= num / 2; ++i) {
    lengths[iter] = i;
    if (i == 1 && iter >= n / 3) {
      continue;
    }
    // if (i == 1) {
    //   genTree(iter + 1, num - i, 2);
    // } else {
    genTree(iter + 1, num - i, i);
    // }
  }
  lengths[iter] = num;
  genTree(iter + 1, 0, 0);
}

void spiderTree(int num) {
  genTree(0, num - 1, 1);
}

int main(int argc, char* argv[]) {
  beginClock = clock();
  prevClock = beginClock;
  srand(time(NULL));

  n = stoi(argv[1]);
  if (n % 2 == 0) {
    cerr << "NOPE" << endl;
    return 0;
  }
  spiderTree(n); // without paths
  cerr << "total amount: " << totalAmount << endl;

  int prevn = 0;
  int tree_count = 0;
  for (int ttt = 0; ttt < totalAmount; ttt++) {
    int tt = ttt;

    if (n > prevn) {
      cout << "n=" << n << endl;
      prevn = n;
    }

    for (int v = 0; v < n; ++v) {
      parent[v] = treeDatabase[tt][v + 1];
      curTreeLengths[v] = treeLengths[tt][v];
    }

    for (int i = 0; i < n; ++i) {
      deg[i] = 0;
      leaf_count[i] = 0;
    }

    for (int i = 1; i < n; ++i) {
      ++deg[parent[i]];
      ++deg[i];
    }
    assert(deg[0] > 1);

    /*bool is_reducible = false;
    for (int i = 1; i < n; ++i) {
      if (deg[i] == 1) {
        ++leaf_count[parent[i]];
        if (leaf_count[parent[i]] > 1) {
          is_reducible = true;
          break;
        }
      }
    }
    if (is_reducible) {
      continue;
    }
    bool is_even_spider = true;
    int max_deg = 0;
    for (int i = 0; i < n; ++i) {
      if (deg[i] > 2) {
        if (max_deg != 0) {
          is_even_spider = false;
          break;
        }
        max_deg = deg[i];
        if (max_deg % 2 != 0) {
          is_even_spider = false;
          break;
        }
      }
    }
    if (!is_even_spider) {
      continue;
    }
    // filter paths
    if (max_deg == 0) {
      continue;
    }*/
    ++tree_count;

    if (tree_count % 10 == 0) {
      cerr << tree_count << ",";
    }

    int order_pos = 0;
    for (int p = 0; p < n; ++p) {
      for (int v = 1; v < n; ++v) {
        if (parent[v] == p) {
          order[order_pos] = v;
          ++order_pos;
        }
      }
    }

    for (int v = 0; v < n; ++v) {
      children[v].clear();
    }
    for (int v = 1; v < n; ++v) {
      children[parent[v]].push_back(v);
    }

    /*type[0] = 0;
    for (int i = 1; i < n; ++i)
      type[i] = 1 - type[parent[i]];*/

    for (int i = 0; i < n; ++i) {
      gr[i].clear();
    }

    for (int i = 1; i < n; ++i) {
      gr[parent[i]].push_back(i);
      gr[i].push_back(parent[i]);
    }

    edgeGracefulLabels.clear();
    buildLinks();

    for (int i = 0; i < n; ++i) {
      usedVal[i] = false;
      cur_sum[i] = 0;

      count_used[i] = 0;
      for (int j = 0; j < n; ++j) {
        used_e_at_v[i][j] = false;
      }
    }
    for (int i = 1; i <= n / 2; ++i) {
      edges[2 * i - 2] = i;
      edges[2 * i - 1] = -i;
    }

    sol_count = 0;
    all_profiles.clear();
    solved_profiles.clear();
    prof_count.clear();
    sols.clear();
    profiles.clear();
    schemes.clear();
    scheme_count.clear();
    all_lengths_different = true;
    for (int i = 1; i < deg[0]; ++i) {
      if (curTreeLengths[order[i - 1]] == curTreeLengths[order[i]]) {
        all_lengths_different = false;
        break;
      }
    }
    // if (!all_lengths_different) {
    //   continue;
    // }

    // if (curTreeLengths[order[deg[0] - 1]] != (n - 1) / 2 - 1) {
    //   continue;
    // }
    res = gen(1);
    if (all_profiles.empty()) {
      continue;
    }

    if (all_profiles.size() != solved_profiles.size()) {
      cout << "something bad happened" << endl;
      cerr << "something bad happened" << endl;
      for (const auto& prof : all_profiles) {
        if (solved_profiles.find(prof) == solved_profiles.end()) {
          cout << "missed profile: " << prof << endl;
          cerr << "missed profile: " << prof << endl;
        }
      }
    }

    // print out all labelings
    for (const auto& prof : solved_profiles) {
      cout << "|" << endl;
      for (int i = 0; i < edgeGracefulLabels.size(); ++i) {
        if (profiles[i] != prof) {
          continue;
        }
        cout << "edgeGraceful labels: ";
        for (int j = 1; j < n; ++j) {
          cout << edgeGracefulLabels[i][j - 1];
          if (j < n - 1 && parent[j + 1] == 0) {
            cout << "; ";
          }
          cout << " ";
        }
        cout << endl;

        cout << "| ";
        if (!all_lengths_different) {
          cout << "prof: " << profiles[i] << " ";
        }
        if (all_lengths_different) {
          cout << "prof: " << schemes[i];
        }
        cout << "vertices: ";
        for (int j = 1; j < n; ++j) {
          if (deg[j] == 1) {
            cout << edgeGracefulLabels[i][j - 1] << "; ";
          } else {
            cout << edgeGracefulLabels[i][j - 1] + edgeGracefulLabels[i][j] << " ";
          }
        }
        cout << endl;
      }
    }

    cout << "tree: ";
    for (int j = 1; j < n; ++j) {
      cout << parent[j] << "->" << j << "; ";
    }
    cout << endl;
    countDiam();
    int maxDeg = 0;
    for (int i = 0; i < n; ++i) {
      maxDeg = max(maxDeg, deg[i]);
    }
    cout << "order: ";
    for (int i = 0; i < n - 1; ++i) {
      cout << order[i] << " ";
    }
    cout << endl;
    cout << "edges: ";
    for (int i = 0; i < n - 1; ++i) {
      cout << edges[i] << " ";
    }
    cout << endl;
    cout << "diam=" << diam << "; maxDeg=" << maxDeg << endl;
    cout << "sol_count=" << sol_count << ";" << endl;
    cout << "cnts: " << endl;
    for (const auto& prof : prof_count) {
      cout << prof.first << "  cnt: " << prof.second << endl;
    }
    if (all_lengths_different) {
      cout << "scheme counts" << endl;
      for (const auto& scheme : scheme_count) {
        cout << scheme.first << "  count: " << scheme.second << endl;
      }
    }
    cout << endl;
    cout << endl;
  }

  endClock = clock();
  double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;

  cerr << "Time: " << elapsed_secs  << "s" << endl;
  cout << "Time: " << elapsed_secs  << "s" << endl;
  cerr << "tree count: " << tree_count << endl;
  cerr << "the end" << endl;
  //system ("pause");
  return 0;
}
