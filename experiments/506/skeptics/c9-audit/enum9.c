/* enum9.c -- from-scratch LABELLED enumeration of abstract Moebius block structures on 9 points
 * (hostile-referee audit of the c(9) = 25 claim; independent of the theory agent's code).
 *
 * Model.  P = {0..8}.  A "block" is a subset of size >= 4 (a circle-or-line with >= 4 points);
 * a triple not inside any block is a 3-block.  Necessary conditions used:
 *   (I)  two blocks share <= 2 points;
 *   (II) pair coverage at a point p by blocks of size >= 4 through p is <= cap_point
 *        (28 = C(8,2) is the trivial combinatorial cap; 27 = Sylvester-Gallai only; 24 = o(8)=4);
 *   (III) lines = blocks-or-3-blocks pairwise sharing <= 1 point, total pair coverage <= cap_lines
 *        (36 trivial; 35 SG only; 30 with o(9)=6).
 * circles = C(9,3) - D - ell,  D = sum over blocks (C(k,3)-1).
 *
 * Symmetry breaking: the LARGEST block has size m and is relabelled to {0,...,m-1}; all other blocks
 * have size 4..m.  Every isomorphism class of structures with largest block size m is visited (as
 * one or more labelled copies).  No canonicity test is used, so no orderly-generation lemma is needed.
 *
 * Usage: enum9 m cap_point cap_lines target [maxnodes]
 * Output: one line per labelled family with C(9,3)-D-ellmax <= target:
 *   FAM D=.. ell=.. count=.. degs=[..] blocks=[[..],[..],...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef N
#define N 9
#endif
#define NT (N*(N-1)*(N-2)/6)
#define MAXC 512

static int cand[MAXC], ncand;          /* candidate block codes (sorted increasing), forced block excluded */
static int cdef[MAXC], cinc[MAXC];     /* deficit C(k,3)-1 and per-point coverage increment C(k-1,2) */
static int F[64], nF;                  /* current family */
static int cov[N];
static int triples[NT];
static int cap_point, cap_lines, target, Dmin, ellcap, m;
static long long nodes = 0, leaves = 0, found = 0, maxnodes = 0;
static int best;

static int binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    long r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return (int)r;
}
static int popc(int x) { return __builtin_popcount(x); }

/* ---- maximum number of lines for the current family ---- */
static int lc[256], lp[256], nl;
static void lrec(const int *avail, int na, int size, int used) {
    if (size > best) best = size;
    if (size + na <= best) return;
    if (size + (cap_lines - used) / 3 <= best) return;
    int nb[256];
    for (int i = 0; i < na; i++) {
        int j = avail[i];
        if (size + (na - i) <= best) return;
        if (used + lp[j] > cap_lines) continue;
        int nnb = 0;
        for (int k = i + 1; k < na; k++) {
            int q = avail[k];
            if (popc(lc[j] & lc[q]) <= 1) nb[nnb++] = q;
        }
        lrec(nb, nnb, size + 1, used + lp[j]);
    }
}
static int ellmax(void) {
    nl = 0;
    for (int i = 0; i < nF; i++) { lc[nl] = F[i]; int k = popc(F[i]); lp[nl] = k * (k - 1) / 2; nl++; }
    for (int t = 0; t < NT; t++) {
        int c = triples[t], covered = 0;
        for (int i = 0; i < nF; i++) if ((F[i] & c) == c) { covered = 1; break; }
        if (!covered) { lc[nl] = c; lp[nl] = 3; nl++; }
    }
    int avail[256];
    for (int i = 0; i < nl; i++) avail[i] = i;
    best = 0;
    lrec(avail, nl, 0, 0);
    return best;
}

static void print_family(int D, int ell) {
    printf("FAM D=%d ell=%d count=%d degs=[", D, ell, NT - D - ell);
    for (int p = 0; p < N; p++) {
        int d = 0;
        for (int i = 0; i < nF; i++) if (F[i] >> p & 1) d++;
        printf("%d%s", d, p < N - 1 ? "," : "");
    }
    printf("] blocks=[");
    for (int i = 0; i < nF; i++) {
        printf("[");
        int first = 1;
        for (int p = 0; p < N; p++) if (F[i] >> p & 1) { printf("%s%d", first ? "" : ",", p); first = 0; }
        printf("]%s", i < nF - 1 ? "," : "");
    }
    printf("]\n");
}

static void dfs(int start, int D) {
    nodes++;
    if (maxnodes && nodes > maxnodes) return;
    if (D >= Dmin) {
        leaves++;
        int l = ellmax();
        if (NT - D - l <= target) { found++; print_family(D, l); fflush(stdout); }
    }
    int ext[MAXC], ne = 0;
    long sumdef = 0;
    for (int i = start; i < ncand; i++) {
        int b = cand[i], ok = 1;
        for (int j = 0; j < nF; j++) if (popc(F[j] & b) > 2) { ok = 0; break; }
        if (!ok) continue;
        for (int p = 0; p < N; p++) if ((b >> p & 1) && cov[p] + cinc[i] > cap_point) { ok = 0; break; }
        if (!ok) continue;
        ext[ne++] = i; sumdef += cdef[i];
    }
    if (D + sumdef < Dmin) return;
    int rem = 0;
    for (int p = 0; p < N; p++) rem += cap_point - cov[p];
    if (D + rem / 3 < Dmin) return;     /* every block of size k<=8 has deficit <= (1/3) * (its coverage units) */
    long suffix = sumdef;
    for (int e = 0; e < ne; e++) {
        int i = ext[e];
        if (D + suffix < Dmin) break;
        int b = cand[i];
        F[nF++] = b;
        for (int p = 0; p < N; p++) if (b >> p & 1) cov[p] += cinc[i];
        dfs(i + 1, D + cdef[i]);
        for (int p = 0; p < N; p++) if (b >> p & 1) cov[p] -= cinc[i];
        nF--;
        suffix -= cdef[i];
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: enum9 m cap_point cap_lines target [maxnodes]\n"); return 1; }
    m = atoi(argv[1]); cap_point = atoi(argv[2]); cap_lines = atoi(argv[3]); target = atoi(argv[4]);
    if (argc > 5) maxnodes = atoll(argv[5]);
    ellcap = cap_lines / 3;
    Dmin = NT - target - ellcap;
    int nt = 0;
    for (int c = 0; c < (1 << N); c++) if (popc(c) == 3) triples[nt++] = c;
    int forced = (1 << m) - 1;
    ncand = 0;
    for (int c = 0; c < (1 << N); c++) {
        int k = popc(c);
        if (k >= 4 && k <= m && c != forced) {
            cand[ncand] = c; cdef[ncand] = binom(k, 3) - 1; cinc[ncand] = binom(k - 1, 2); ncand++;
        }
    }
    nF = 0; memset(cov, 0, sizeof cov);
    F[nF++] = forced;
    for (int p = 0; p < m; p++) cov[p] += binom(m - 1, 2);
    fprintf(stderr, "enum9: m=%d forced={0..%d} cap_point=%d cap_lines=%d target=%d ellcap=%d Dmin=%d ncand=%d\n",
            m, m - 1, cap_point, cap_lines, target, ellcap, Dmin, ncand);
    clock_t t0 = clock();
    dfs(0, binom(m, 3) - 1);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "enum9: done m=%d nodes=%lld leaves(D>=Dmin)=%lld found=%lld time=%.1fs%s\n",
            m, nodes, leaves, found, secs, (maxnodes && nodes > maxnodes) ? "  *** NODE LIMIT HIT: INCOMPLETE ***" : "");
    printf("END m=%d nodes=%lld leaves=%lld found=%lld complete=%d\n", m, nodes, leaves, found,
           !(maxnodes && nodes > maxnodes));
    return 0;
}
