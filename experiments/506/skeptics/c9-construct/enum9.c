/* Independent enumeration of abstract Mobius block structures on n=9 points (c9 audit).
   A "structure" is a family F of blocks (subsets of size 4..8) pairwise sharing <= 2 points.
   The largest block is fixed WLOG to {0,..,s-1} (s given on the command line); all other blocks
   have size in [4, s] and meet the fixed block in <= 2 points.  Constraints used (all provably
   necessary for a real point set, using only Sylvester-Gallai o(m) >= 1):
     (C1) blocks pairwise share <= 2 points;
     (C2) for every point p: sum_{B ∋ p} C(|B|-1,2) <= C(8,2) - 1 = 27;
     (C3) lines: a set L of blocks/uncovered triples pairwise sharing <= 1 point with
          sum_{L} C(|L|,2) <= C(9,2) - 1 = 35.
   For every family with D = sum (C(|B|,3)-1) >= Dmin we compute ell_max = max |L| under (C3) and
   print the family if C(9,3) - D - ell_max <= target.  No symmetry breaking beyond fixing the
   largest block, so every labelled family (with that block as the lexicographically first
   largest block... no: with *some* largest block equal to {0..s-1}) is visited. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 9
#define FULL ((1<<N)-1)
static int pc[1<<N];
static int C3[N+1], C2[N+1];
static int s, target, Dmin;
static int cand[600], ncand;
static int F[64], nF;
static int cov[N];
static long nodes=0, families=0, printed=0;
static int bigblock;

/* ---- ell_max: maximum set of pairwise <=1-sharing "lines" among candidates lc[] ---- */
static int lc[200], nlc, best;
static int chosen[64];
static void lines_dfs(int idx, int cnt, int pairs){
    if(cnt>best) best=cnt;
    if(cnt + (nlc-idx) <= best) return;
    for(int j=idx;j<nlc;j++){
        int L=lc[j]; int ok=1;
        for(int t=0;t<cnt;t++) if(pc[chosen[t]&L]>1){ok=0;break;}
        if(!ok) continue;
        int np=pairs+C2[pc[L]];
        if(np>35) continue;
        chosen[cnt]=L;
        lines_dfs(j+1,cnt+1,np);
    }
}
static int ell_max(void){
    nlc=0;
    for(int i=0;i<nF;i++) lc[nlc++]=F[i];
    for(int T=0;T<=FULL;T++) if(pc[T]==3){
        int covered=0;
        for(int i=0;i<nF;i++) if((F[i]&T)==T){covered=1;break;}
        if(!covered) lc[nlc++]=T;
    }
    best=0;
    lines_dfs(0,0,0);
    return best;
}
static void print_family(int D,int ell){
    printf("FAM D=%d ell=%d count=%d nF=%d blocks:",D,ell,C3[N]-D-ell,nF);
    for(int i=0;i<nF;i++){
        printf(" ");
        for(int p=0;p<N;p++) if(F[i]>>p&1) printf("%d",p);
    }
    printf("\n");
}
static int compatible(int b){
    for(int i=0;i<nF;i++) if(pc[F[i]&b]>2) return 0;
    for(int p=0;p<N;p++) if(b>>p&1) if(cov[p]+C2[pc[b]-1]>27) return 0;
    return 1;
}
static void dfs(int idx,int D){
    nodes++;
    if(D>=Dmin){
        families++;
        int ell=ell_max();
        if(C3[N]-D-ell<=target){ printed++; print_family(D,ell); }
    }
    /* potential bound */
    int remcap=0; for(int p=0;p<N;p++) remcap+=27-cov[p];
    int potdef=0;
    for(int j=idx;j<ncand;j++) if(compatible(cand[j])) potdef+=C3[pc[cand[j]]]-1;
    int pot = potdef; if(remcap/3<pot) pot=remcap/3;
    if(D+pot<Dmin) return;
    for(int j=idx;j<ncand;j++){
        int b=cand[j];
        if(!compatible(b)) continue;
        F[nF++]=b;
        for(int p=0;p<N;p++) if(b>>p&1) cov[p]+=C2[pc[b]-1];
        dfs(j+1,D+C3[pc[b]]-1);
        for(int p=0;p<N;p++) if(b>>p&1) cov[p]-=C2[pc[b]-1];
        nF--;
    }
}
int main(int argc,char**argv){
    s=atoi(argv[1]); target=atoi(argv[2]); int ellcap=atoi(argv[3]); /* ell upper bound used only for Dmin */
    for(int i=0;i<(1<<N);i++) pc[i]=__builtin_popcount(i);
    for(int k=0;k<=N;k++){ C2[k]=k*(k-1)/2; C3[k]=k*(k-1)*(k-2)/6; }
    Dmin = C3[N]-target-ellcap;
    bigblock=(1<<s)-1;
    ncand=0;
    for(int b=0;b<=FULL;b++){
        if(pc[b]<4||pc[b]>s||b==bigblock) continue;
        if(pc[b&bigblock]>2) continue;
        cand[ncand++]=b;
    }
    nF=0; memset(cov,0,sizeof cov);
    F[nF++]=bigblock; for(int p=0;p<s;p++) cov[p]+=C2[s-1];
    fprintf(stderr,"s=%d target=%d ellcap=%d Dmin=%d ncand=%d\n",s,target,ellcap,Dmin,ncand);
    dfs(0,C3[s]-1);
    fprintf(stderr,"nodes=%ld families(D>=Dmin)=%ld printed=%ld\n",nodes,families,printed);
    return 0;
}
