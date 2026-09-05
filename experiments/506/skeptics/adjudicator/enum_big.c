/* Adjudicator's own enumerator of abstract Moebius block structures (written from scratch).
   usage: enum_big n m target capmode
   Points 0..n-1.  The largest block is forced to be {0..m-1}; every other block has size 4..m,
   meets the forced block in <= 2 points, and blocks pairwise share <= 2 points (Lemma 2.1).
   capmode 0: no further constraints.  capmode 1: Sylvester-Gallai caps with o = 1:
       for every point p:  sum_{B ∋ p} C(|B|-1,2) <= C(n-1,2) - 1   (derived structure at p has an
       ordinary line);  lines cover <= C(n,2) - 1 pairs.
   Plain labelled DFS (blocks added in increasing candidate index; every compatible family is a
   node).  At every node with D >= Dmin = C(n,3) - target - ellcap we compute ell_max exactly
   (maximum set of pairwise <=1-intersecting "lines" chosen among the family's blocks and the
   triples not inside any block, within the pair budget) and print the family if
   C(n,3) - D - ell_max <= target.  Pruning only discards subtrees whose every descendant has
   D < Dmin (two upper bounds on the additional deficit), so the output is exhaustive. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, m, target, capmode;
static int C3n, cap_point, cap_lines, ellcap, Dmin;
static int binom[24][24];
#define MAXC 1024
static int ncand; static uint32_t cand[MAXC]; static int cdef[MAXC], ccov[MAXC];
static uint8_t compat[MAXC][MAXC];
static uint32_t fam[64]; static int fidx[64]; static int nf; static int D; static int cov[24];
static long long nodes = 0, checked = 0, found = 0;
static int ratio_num, ratio_den;

static int popc(uint32_t x) { return __builtin_popcount(x); }

/* ---- maximum line set (clique in the <=1-intersection graph, with pair budget) ---- */
#define MAXL 1200
static uint32_t lc[MAXL]; static int lpairs[MAXL]; static int nl;
static uint8_t lcomp[MAXL][MAXL];
static int best_ell; static int sel[64];

static void clique(int start, int size, int pairs, int nsel) {
    if (size > best_ell) best_ell = size;
    for (int i = start; i < nl; i++) {
        if (size + (nl - i) <= best_ell) return;
        if (size + (cap_lines - pairs) / 3 <= best_ell) return;
        if (pairs + lpairs[i] > cap_lines) continue;
        int ok = 1;
        for (int j = 0; j < nsel; j++) if (!lcomp[sel[j]][i]) { ok = 0; break; }
        if (!ok) continue;
        sel[nsel] = i;
        clique(i + 1, size + 1, pairs + lpairs[i], nsel + 1);
    }
}

static int ell_max(void) {
    nl = 0;
    for (int i = 0; i < nf; i++) { lc[nl] = fam[i]; lpairs[nl] = binom[popc(fam[i])][2]; nl++; }
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) {
        uint32_t t = (1u << a) | (1u << b) | (1u << c);
        int inside = 0;
        for (int i = 0; i < nf; i++) if ((fam[i] & t) == t) { inside = 1; break; }
        if (!inside) { lc[nl] = t; lpairs[nl] = 3; nl++; }
    }
    for (int i = 0; i < nl; i++) for (int j = 0; j < nl; j++) lcomp[i][j] = (popc(lc[i] & lc[j]) <= 1);
    best_ell = 0;
    clique(0, 0, 0, 0);
    return best_ell;
}

static void print_family(int ell) {
    printf("FOUND D=%d ell=%d count=%d blocks=", D, ell, C3n - D - ell);
    for (int i = 0; i < nf; i++) {
        printf("[");
        int first = 1;
        for (int p = 0; p < n; p++) if (fam[i] >> p & 1) { printf(first ? "%d" : ",%d", p); first = 0; }
        printf("]");
    }
    printf("\n");
    fflush(stdout);
}

static void dfs(int start) {
    nodes++;
    if (D >= Dmin) {
        checked++;
        int e = ell_max();
        if (C3n - D - e <= target) { found++; print_family(e); }
    }
    int rem = 0;
    for (int p = 0; p < n; p++) rem += cap_point - cov[p];
    long long pot1 = (long long)rem * ratio_num / ratio_den;
    long long pot2 = 0;
    for (int i = start; i < ncand; i++) {
        int ok = 1;
        for (int j = 0; j < nf; j++) if (!compat[fidx[j]][i]) { ok = 0; break; }
        if (ok) pot2 += cdef[i];
    }
    long long pot = pot1 < pot2 ? pot1 : pot2;
    if (D + pot < Dmin) return;
    for (int i = start; i < ncand; i++) {
        int ok = 1;
        for (int j = 0; j < nf; j++) if (!compat[fidx[j]][i]) { ok = 0; break; }
        if (!ok) continue;
        uint32_t B = cand[i]; int cp = ccov[i];
        for (int p = 0; p < n; p++) if ((B >> p & 1) && cov[p] + cp > cap_point) { ok = 0; break; }
        if (!ok) continue;
        fam[nf] = B; fidx[nf] = i; nf++; D += cdef[i];
        for (int p = 0; p < n; p++) if (B >> p & 1) cov[p] += cp;
        dfs(i + 1);
        for (int p = 0; p < n; p++) if (B >> p & 1) cov[p] -= cp;
        nf--; D -= cdef[i];
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: enum_big n m target capmode\n"); return 1; }
    n = atoi(argv[1]); m = atoi(argv[2]); target = atoi(argv[3]); capmode = atoi(argv[4]);
    for (int i = 0; i < 24; i++) { binom[i][0] = 1; for (int j = 1; j <= i; j++) binom[i][j] = binom[i-1][j-1] + (j <= i-1 ? binom[i-1][j] : 0); }
    C3n = binom[n][3];
    cap_point = binom[n-1][2] - (capmode == 1 ? 1 : 0);
    cap_lines = binom[n][2] - (capmode == 1 ? 1 : 0);
    ellcap = cap_lines / 3;
    Dmin = C3n - target - ellcap;
    ratio_num = binom[m][3] - 1; ratio_den = m * binom[m-1][2];
    uint32_t forced = (1u << m) - 1;
    ncand = 0;
    cand[ncand] = forced; cdef[ncand] = binom[m][3] - 1; ccov[ncand] = binom[m-1][2]; ncand++;
    for (uint32_t S = 1; S < (1u << n); S++) {
        int k = popc(S);
        if (k < 4 || k > m || S == forced) continue;
        if (popc(S & forced) > 2) continue;
        cand[ncand] = S; cdef[ncand] = binom[k][3] - 1; ccov[ncand] = binom[k-1][2]; ncand++;
    }
    for (int i = 0; i < ncand; i++) for (int j = 0; j < ncand; j++) compat[i][j] = (popc(cand[i] & cand[j]) <= 2);
    fprintf(stderr, "enum_big n=%d m=%d target=%d capmode=%d: C3n=%d cap_point=%d cap_lines=%d ellcap=%d Dmin=%d ncand=%d ratio=%d/%d\n",
            n, m, target, capmode, C3n, cap_point, cap_lines, ellcap, Dmin, ncand, ratio_num, ratio_den);
    nf = 1; fam[0] = forced; fidx[0] = 0; D = cdef[0];
    for (int p = 0; p < n; p++) cov[p] = (forced >> p & 1) ? ccov[0] : 0;
    dfs(1);
    printf("END n=%d m=%d target=%d capmode=%d nodes=%lld checked=%lld found=%lld\n", n, m, target, capmode, nodes, checked, found);
    fprintf(stderr, "END n=%d m=%d target=%d capmode=%d nodes=%lld checked=%lld found=%lld\n", n, m, target, capmode, nodes, checked, found);
    return 0;
}
