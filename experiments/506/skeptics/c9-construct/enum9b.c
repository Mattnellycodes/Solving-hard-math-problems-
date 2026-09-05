/* enum9b: same as enum9.c (independent enumeration of abstract Mobius block structures on 9 points,
   largest block fixed WLOG to {0..s-1}, SG-only caps) PLUS the constraint that every point lies in at most
   7 four-blocks.  Validity: the 4-blocks through p give 3-point lines of the 8-point derived set (after
   inversion about p), pairwise sharing <=1 point; 8 such lines would form the Mobius-Kantor (8_3)
   configuration, which has no real realisation (checked exactly in mk83.py).  A degree-based potential
   bound is added so that the s=4 case (17 four-blocks needed) is pruned immediately. */
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
static int cov[N], deg4[N];
static long nodes=0, families=0, printed=0;
static int bigblock;
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
    best=0; lines_dfs(0,0,0); return best;
}
static void print_family(int D,int ell){
    printf("FAM D=%d ell=%d count=%d nF=%d blocks:",D,ell,C3[N]-D-ell,nF);
    for(int i=0;i<nF;i++){ printf(" "); for(int p=0;p<N;p++) if(F[i]>>p&1) printf("%d",p); }
    printf("\n"); fflush(stdout);
}
static int compatible(int b){
    for(int i=0;i<nF;i++) if(pc[F[i]&b]>2) return 0;
    for(int p=0;p<N;p++) if(b>>p&1){
        if(cov[p]+C2[pc[b]-1]>27) return 0;
        if(pc[b]==4 && deg4[p]>=7) return 0;
    }
    return 1;
}
static void dfs(int idx,int D){
    nodes++;
    if(D>=Dmin){
        families++;
        int ell=ell_max();
        if(C3[N]-D-ell<=target){ printed++; print_family(D,ell); }
    }
    int remcap=0, remdeg=0;
    for(int p=0;p<N;p++){ remcap+=27-cov[p]; remdeg+=7-deg4[p]; }
    int potdef=0, pot4=0, potbig=0;
    for(int j=idx;j<ncand;j++) if(compatible(cand[j])){
        int k=pc[cand[j]];
        if(k==4){ pot4++; } else potbig+=C3[k]-1;
    }
    /* deficit from 4-blocks is at most 3*min(pot4, remdeg/4); big blocks bounded by their list */
    int p4 = pot4; if(remdeg/4<p4) p4=remdeg/4;
    potdef = 3*p4 + potbig;
    int pot = potdef; if(remcap/3<pot) pot=remcap/3;
    if(D+pot<Dmin) return;
    for(int j=idx;j<ncand;j++){
        int b=cand[j];
        if(!compatible(b)) continue;
        F[nF++]=b;
        for(int p=0;p<N;p++) if(b>>p&1){ cov[p]+=C2[pc[b]-1]; if(pc[b]==4) deg4[p]++; }
        dfs(j+1,D+C3[pc[b]]-1);
        for(int p=0;p<N;p++) if(b>>p&1){ cov[p]-=C2[pc[b]-1]; if(pc[b]==4) deg4[p]--; }
        nF--;
    }
}
int main(int argc,char**argv){
    s=atoi(argv[1]); target=atoi(argv[2]); int ellcap=atoi(argv[3]);
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
    nF=0; memset(cov,0,sizeof cov); memset(deg4,0,sizeof deg4);
    F[nF++]=bigblock; for(int p=0;p<s;p++){ cov[p]+=C2[s-1]; if(s==4) deg4[p]++; }
    fprintf(stderr,"s=%d target=%d ellcap=%d Dmin=%d ncand=%d\n",s,target,ellcap,Dmin,ncand);
    dfs(0,C3[s]-1);
    fprintf(stderr,"FINISHED nodes=%ld families(D>=Dmin)=%ld printed=%ld\n",nodes,families,printed);
    return 0;
}
