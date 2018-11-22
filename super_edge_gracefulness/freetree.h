/*==============================================================*/
/* program: freetree.h                                          */
/* purpose: generating all free trees                           */
/* input  : n -- number of nodes                                */
/*          m -- max degree                                     */
/*          lb, ub -- lower and upper bound on diameter         */
/* output : listing of free trees in relex order                */
/* date   : September 1995, updated September 2000              */
/* programmers: Gang Li & Frank Ruskey                          */
/* algorithm: From the paper: G. Li and F. Ruskey, "The         */
/*    Advantages of Forward Thinking in Generating Rooted and   */
/*    Free Trees",  10th Annual ACM-SIAM Symposium on Discrete  */
/*    Algorithms (SODA), (1999) S939-940.  See the web page at  */
/*    http://webhome.cs.uvic.ca/~ruskey/fruskey.html            */
/*      Publications/RootedFreeTree.html                        */
/* more info: see                                               */
/*    http://theory.cs.uvic.ca/inf/tree/FreeTrees.html          */
/*==============================================================*/

// TODO: Put into namespace
const int MAX_SIZE = 50;      /* max size of the tree           */

int par[MAX_SIZE];            /* parent position of i           */
int L[MAX_SIZE];              /* level of node i                */
int k;                        /* max number of children         */
int chi[MAX_SIZE];            /* number of children of a node   */
int nextp[MAX_SIZE];          /* next good pos to add nodes     */
int rChi[MAX_SIZE];           /* the right most child of node i */
int ub;                       /* upper bound                    */
int num;                      /* number of trees                */

int outputP, outputL;

int N;  /* Corresponds to -n */
int K;  /*                -k */
int M;            /*      -m */
int out_format;   /*      -o */
int U;            /*      -u */

short treeDatabase[30000000][30];
static int totalAmount = 0;

void PrintIt() {
    ++num;
    if (outputP) {
        for (int i = 1; i <= N; ++i) {
            printf("%d", par[i]);
            if (N > 9 && i < N) {
                printf(",");
            }
        }
    }
    if (outputL) {
        if (outputL && outputP)
            printf(" : ");
        for (int i = 1; i <= N; ++i) {
            printf("%d", L[i]);
            if (N > 9 && i < N) {
                printf(",");
            }
        }
    }
    if (outputP || outputL) {
        printf("\n");
    }
}

void updateL() {
    L[1] = 0;
    for (int i = 2; i <= N; ++i) {
        L[i] = L[par[i]] + 1;
    }
}

bool good(int p, int h, int t) {
    if (p == 2 && K <= 2 && t == 0)
        return true;
    if (t == 1) {
        if ((2 * h >= K + 1) && (2 * h <= U + 1)) {
            if ((p - 1) * 2 >= N) {
                return true;
            } else if (p - h - 1 == 1) {
                if (par[p] > 2) {
                    return true;
                }
            } else if ((p - h - 1 >= 2) && ((par[h + 2] > 2) || (par[h + 3] > 2))) {
                return true;
            }
        }
    } else if ((N - p >= h) && (2 * h >= K)) {
        if ((U == N - 1) && (N % 2 == 0)) {
            return (2 * h <= U + 1);
        } else {
            return(2 * h <= U);
        }
    }
    return false;
}

void GenTree(int p, int s, int cL, int h, int l, int n, int f, int g) {
    int hh, flag, entry, temp;
    if (p > n)
        if (f == 0) {
            if (good(p - 1, h, 0)) {
                GenTree(p, 2, p - 2, h, n, N, 1, 0);
            }
            if (good(p - 1, h, 1)) {
                GenTree(p, 2, p - 3, h, n, N, 1, 1);
            }
        }
        else {
            updateL();
            //PrintIt();
            treeDatabase[totalAmount][0] = N;
            for (size_t i = 1; i <= static_cast<size_t>(N); ++i)
                treeDatabase[totalAmount][i] = par[i] - 1;
            ++totalAmount;
        }
    else {
        if (cL == 0) {
            if (p < ub + 2) {
                par[p] = p - 1;
            } else {
                GenTree(p, p - 1, 1, h, l, n, f, g);
                return;
            }
        } else if (par[p - cL] < s) {
            par[p] = par[p - cL];
        } else {
            par[p] = cL + par[p - cL];
            if (g == 1) {
                if (((l - 1) * 2 < n) && (p - cL <= l) && (
                    ((p - cL + 1 < l) &&  (par[p - cL + 1] == 2)
                    && (p - cL + 2 <= l) && (par[p - cL + 2] == 2))     /* case 1 */
                    || ((p - cL + 1 == l) && (par[p - cL + 1] == 2))    /* case 2 */
                    || (p - cL + 1 > l)))                               /* case 3 */
                {
                    s = par[p];
                    cL = p - s;
                    par[p] = par[par[p]];
                } else {
                    if (par[p - cL] == 2) {
                        par[p] = 1;
                    }
                }
            }
        }
        if (s != 0 || p <= ub + 1) {
            chi[par[p]] = chi[par[p]] + 1;
            temp = rChi[par[p]];
            rChi[par[p]] = p;
            if (chi[par[p]] <= ((par[p] == 1) ? k : k - 1)) {
                if (chi[par[p]] < (par[p] == 1 ? k : k - 1))
                    nextp[p] = par[p];
                else
                    nextp[p] = nextp[par[p]];
                GenTree(p + 1, s, cL, h, l, n, f, g);
            }
            chi[par[p]] = chi[par[p]] - 1;
            rChi[par[p]] = temp;
        }
        if (s == 0 && 2 * (p - 2) < K)
            return;

        nextp[p] = nextp[par[p]];
        entry = nextp[p];
        flag = 0;
        hh = 1;
        while ((((f == 0) && (entry >= 2)) || ((f == 1) && (entry >= 1))) && (flag == 0)) {
            if (s == 0)
                h = p - 2;
            if (p <= l + h - g)
                hh = 0;
            if ((f == 0) || (hh == 1)) {
                //s = par[p];
                //par[p] = par[s];
                par[p] = entry;

                chi[entry] = chi[entry] + 1;
                temp = rChi[par[p]];
                rChi[par[p]] = p;
                if (chi[entry] >= (entry == 1 ? k : k - 1))
                    nextp[p] = nextp[entry];
                if (f == 0)
                    GenTree(p + 1, temp, p - temp, h, 0, N - h + 1, f, g);
                else if (hh == 1)
                    GenTree(p + 1, temp, p - temp, h, l, n, f, g);
                chi[entry] = chi[entry] - 1;
                rChi[par[p]] = temp;
                entry = nextp[entry];
                nextp[p] = entry;
            } else {
                flag = 1;
            }
        }
        if (f == 0) {
            if (good(p - 1, h, 0)) {
                GenTree(p, 2, p - 2, h, p - 1, N, 1, 0);
            }
            if (good(p - 1, h, 1)) {
                GenTree(p, 2, p - 3, h, p - 1, N, 1, 1);
            }
        }
    }
} /* GenTree */

void ProcessInput(int n) {
    /*
    printf("Usage: freetree -n N [-d D] [-l L] [-u U] [-o O]\n");
    printf("  where N is the number of nodes in the tree,\n");
    printf("  D is an upper bound on the degree of any node,\n");
    printf("  L is a lower bound on the diameter of the tree,\n");
    printf("  U is an upper bound on the diameter of the tree,\n");
    printf("  O specifies output:\n");
    printf("     O = 0 output the number of trees only.\n");
    printf("     O = 1 output parent array,\n");
    printf("     O = 2 output level array,\n");
    printf("     O = 3 output both,\n");
    printf("  Default values: K = U = N - 1, L = 2, O = 3.\n");
    printf("  Frank Ruskey, 1995, 2000\n.");
    */
    N = n;
    M = N - 1;
    K = 2;
    U = 4;
    U = N - 1;
    out_format = 1;
    out_format = 0;
}

int freeTree(int n) {
    num = 0;
    /* first set all the parameters */
    ProcessInput(n);
    if (out_format & 1)
        outputP = 1;
    if (out_format & 2)
        outputL = 1;
    if (K > U) {
        //printf("ERROR: lb > ub !!\n");
        return -1;
    }
    if (U > 1) {
        if (U < N) {
            ub = (U + 1) / 2;
        } else {
            ub = N / 2;
            U = N - 1;
        }
    } else {
        //printf("ERROR: ub is too small, should be >= 2\n");
        return -1;
    }

    if (M > 0)
        k = M;
    else
        k = -1;

    if (K <= 1) {
        //printf("ERROR: k is too small, should be >= 2\n");
        return -1;
    } else {
        for (int i = 1; i <= N; i++)
            chi[i] = 0;

        par[1] = 0;
        par[2] = 1;
        nextp[1] = 0;
        nextp[2] = 1;
        chi[1] = 1;
        GenTree(3, 0, 0, ub, 0, N - ub + 1, 0, 0);
        //printf("The total number of trees is : %4d\n", num);
    }
    return 0;
}

