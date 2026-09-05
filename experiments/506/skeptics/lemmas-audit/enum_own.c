/* lemmas-audit: independent labelled enumerator of abstract Moebius block structures (Erdos #506).
 *
 * usage: enum_own n target capmode fixm [maxsize]
 *   capmode: 0 = no Sylvester-Gallai caps (packing constraints only)
 *            1 = weak caps  (o(m) = 1 for all m: pure Sylvester-Gallai)
 *            2 = strong caps (o(3..10) = 3,3,4,3,3,4,6,5)
 *   fixm   : >0: the block {0..fixm-1} is fixed and is the LARGEST block (all other blocks have size
 *            4..fixm and share <= 2 points with it).  0: nothing fixed, all sizes 4..n-1.
 *   maxsize: optional cap on block sizes when fixm = 0 (default n-1).
 * The DFS adds blocks in increasing bitmask order, so every labelled family is visited exactly once.
 * For every family with D >= Dmin it computes ell_max (max number of pairwise <=1-sharing 'lines'
 * chosen among the family's blocks and the uncovered triples, with sum C(|l|,2) <= cap_lines) and
 * prints the family if C(n,3) - D - ell_max <= target.  Output line format:
 *   REC D ell count : mask mask ...
 * Written from scratch; does not share code with the theory agent's enumerators.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, target, capmode, fixm, maxsize;
static int N3, capd, capl, ellmax_global, Dmin;
static int ncand; static int cand[600]; static int candsz[600]; static int canddef[600];
static int F[64]; static int nF; static int cov[16]; static int D;
static long long nodes = 0, leaves = 0, records = 0;
static int binom[20][20];

static int pc(int x){ return __builtin_popcount(x); }
static int olower(int m){
  if (m <= 2) return 0;
  if (capmode == 0) return 0;
  if (capmode == 1) return 1;
  static const int tab[11] = {0,0,0,3,3,4,3,3,4,6,5};
  if (m <= 10) return tab[m];
  return (6*m + 12) / 13; /* ceil(6m/13) */
}

/* ---- ell_max: branch and bound over line candidates ---- */
static int lc[200], lcsz[200], nlc; static int bestell; static int needell;
static void rec_lines(int start, int chosen, int *sel, int pairs){
  if (chosen > bestell) bestell = chosen;
  if (bestell >= needell && needell > 0 && chosen >= needell) return; /* early exit possible */
  /* bound */
  int avail = 0;
  for (int j = start; j < nlc; j++){
    int ok = 1;
    for (int i = 0; i < chosen; i++) if (pc(lc[sel[i]] & lc[j]) > 1){ ok = 0; break; }
    if (ok && pairs + binom[lcsz[j]][2] <= capl) avail++;
  }
  if (chosen + avail <= bestell) return;
  for (int j = start; j < nlc; j++){
    if (pairs + binom[lcsz[j]][2] > capl) continue;
    int ok = 1;
    for (int i = 0; i < chosen; i++) if (pc(lc[sel[i]] & lc[j]) > 1){ ok = 0; break; }
    if (!ok) continue;
    sel[chosen] = j;
    rec_lines(j + 1, chosen + 1, sel, pairs + binom[lcsz[j]][2]);
  }
}
static int ell_max(int need){
  nlc = 0;
  for (int i = 0; i < nF; i++){ lc[nlc] = F[i]; lcsz[nlc] = pc(F[i]); nlc++; }
  /* uncovered triples */
  for (int t = 0; t < (1 << n); t++){
    if (pc(t) != 3) continue;
    int covered = 0;
    for (int i = 0; i < nF; i++) if ((F[i] & t) == t){ covered = 1; break; }
    if (!covered){ lc[nlc] = t; lcsz[nlc] = 3; nlc++; }
  }
  bestell = 0; needell = 0; /* compute the true maximum (no early exit) */
  int sel[200];
  rec_lines(0, 0, sel, 0);
  return bestell;
}

static void process(void){
  leaves++;
  int need = N3 - target - D; if (need < 0) need = 0;
  int e = ell_max(need);
  int count = N3 - D - e;
  if (count <= target){
    records++;
    printf("REC %d %d %d :", D, e, count);
    for (int i = 0; i < nF; i++) printf(" %d", F[i]);
    printf("\n");
    fflush(stdout);
  }
}

static int uncovered_triples(void){
  int u = 0;
  for (int t = 0; t < (1 << n); t++){
    if (pc(t) != 3) continue;
    int covered = 0;
    for (int i = 0; i < nF; i++) if ((F[i] & t) == t){ covered = 1; break; }
    if (!covered) u++;
  }
  return u;
}

static void dfs(int last){
  nodes++;
  if (D >= Dmin) process();
  /* potential bounds */
  long long budget = 0;
  for (int p = 0; p < n; p++) budget += capd - cov[p];
  if (capmode != 0 && D + budget / 3 < Dmin) return;
  /* combinatorial bound: extra deficit <= uncovered triples - (#new blocks) < uncovered triples */
  int U = uncovered_triples();
  if (U == 0) return;
  if (D + U < Dmin) return;
  /* candidate extension list */
  int ext[600]; int next = 0; int potdef = 0;
  for (int c = 0; c < ncand; c++){
    int b = cand[c];
    if (b <= last) continue;
    int ok = 1;
    for (int i = 0; i < nF; i++) if (pc(F[i] & b) > 2){ ok = 0; break; }
    if (!ok) continue;
    int k = candsz[c]; int inc = binom[k-1][2];
    for (int p = 0; p < n; p++) if ((b >> p) & 1) if (cov[p] + inc > capd){ ok = 0; break; }
    if (!ok) continue;
    ext[next++] = c; potdef += canddef[c];
  }
  if (D + potdef < Dmin) return;
  for (int idx = 0; idx < next; idx++){
    /* bound with later candidates only */
    int rest = 0; for (int j = idx; j < next; j++) rest += canddef[ext[j]];
    if (D + rest < Dmin) break;
    int c = ext[idx]; int b = cand[c]; int k = candsz[c]; int inc = binom[k-1][2];
    F[nF++] = b; D += canddef[c];
    for (int p = 0; p < n; p++) if ((b >> p) & 1) cov[p] += inc;
    dfs(b);
    for (int p = 0; p < n; p++) if ((b >> p) & 1) cov[p] -= inc;
    nF--; D -= canddef[c];
  }
}

int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage: enum_own n target capmode fixm [maxsize]\n"); return 1; }
  n = atoi(argv[1]); target = atoi(argv[2]); capmode = atoi(argv[3]); fixm = atoi(argv[4]);
  maxsize = (argc > 5) ? atoi(argv[5]) : n - 1;
  for (int i = 0; i < 20; i++){ binom[i][0] = 1; for (int j = 1; j < 20; j++) binom[i][j] = (i == 0) ? 0 : binom[i-1][j-1] + binom[i-1][j]; }
  N3 = binom[n][3];
  capd = (capmode == 0) ? 1000000 : binom[n-1][2] - olower(n-1);
  capl = (capmode == 0) ? binom[n][2] : binom[n][2] - olower(n);
  ellmax_global = capl / 3;
  Dmin = N3 - target - ellmax_global;
  fprintf(stderr, "n=%d target=%d capmode=%d fixm=%d maxsize=%d: C(n,3)=%d capd=%d capl=%d ellmax<=%d Dmin=%d\n",
          n, target, capmode, fixm, maxsize, N3, capd, capl, ellmax_global, Dmin);
  nF = 0; D = 0; memset(cov, 0, sizeof cov);
  int B0 = 0;
  if (fixm > 0){
    B0 = (1 << fixm) - 1;
    F[nF++] = B0; D += binom[fixm][3] - 1;
    for (int p = 0; p < fixm; p++) cov[p] += binom[fixm-1][2];
    if (capmode != 0 && binom[fixm-1][2] > capd){ fprintf(stderr, "fixed block violates cap\n"); return 0; }
    maxsize = fixm;
  }
  ncand = 0;
  for (int b = 1; b < (1 << n); b++){
    int k = pc(b);
    if (k < 4 || k > maxsize || k > n - 1) continue;
    if (fixm > 0){ if (b == B0) continue; if (pc(b & B0) > 2) continue; }
    cand[ncand] = b; candsz[ncand] = k; canddef[ncand] = binom[k][3] - 1; ncand++;
  }
  fprintf(stderr, "candidates: %d\n", ncand);
  dfs(fixm > 0 ? 0 : 0);
  fprintf(stderr, "nodes=%lld leaves(D>=Dmin)=%lld records=%lld\n", nodes, leaves, records);
  printf("DONE nodes=%lld leaves=%lld records=%lld\n", nodes, leaves, records);
  return 0;
}
