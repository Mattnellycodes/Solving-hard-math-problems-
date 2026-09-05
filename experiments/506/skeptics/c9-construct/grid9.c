/* Exhaustive / random search of 9-point subsets of an integer grid; exact integer circle count.
   usage: grid9 W H mode [samples seed]  mode 0 = exhaustive over all C(W*H,9) subsets, 1 = random samples.
   Prints histogram of circle counts (non-degenerate sets only) and every set with count <= BEST. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a,ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a;}
typedef struct{ll a,b,c,d;} key;
static int W,H,NP; static int px[400],py[400];
static int cmpkey(const void*u,const void*v){const key*x=u,*y=v; if(x->a!=y->a)return x->a<y->a?-1:1; if(x->b!=y->b)return x->b<y->b?-1:1; if(x->c!=y->c)return x->c<y->c?-1:1; if(x->d!=y->d)return x->d<y->d?-1:1; return 0;}
static ll det3(ll m[3][3]){return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);}
static key mk(int i,int j,int k){
    int id[3]={i,j,k}; ll r[3][4];
    for(int t=0;t<3;t++){ll x=px[id[t]],y=py[id[t]]; r[t][0]=x*x+y*y; r[t][1]=x; r[t][2]=y; r[t][3]=1;}
    ll m[3][3]; key K;
    #define SUB(c0,c1,c2) for(int t=0;t<3;t++){m[t][0]=r[t][c0];m[t][1]=r[t][c1];m[t][2]=r[t][c2];}
    SUB(1,2,3) K.a=det3(m); SUB(0,2,3) K.b=-det3(m); SUB(0,1,3) K.c=det3(m); SUB(0,1,2) K.d=-det3(m);
    ll g=gcdll(gcdll(K.a,K.b),gcdll(K.c,K.d)); if(g){K.a/=g;K.b/=g;K.c/=g;K.d/=g;}
    ll lead = K.a?K.a:(K.b?K.b:(K.c?K.c:K.d)); if(lead<0){K.a=-K.a;K.b=-K.b;K.c=-K.c;K.d=-K.d;}
    return K;
}
static int hist[100]; static int bestc=1000;
/* count circles of the 9-set idx[]; returns -1 if degenerate (all on one circle/line) */
static int count9(int*idx,int print){
    key ks[84]; int m=0;
    for(int i=0;i<9;i++)for(int j=i+1;j<9;j++)for(int k=j+1;k<9;k++) ks[m++]=mk(idx[i],idx[j],idx[k]);
    qsort(ks,84,sizeof(key),cmpkey);
    int circles=0,lines=0; int run=1; int maxrun=0;
    for(int t=0;t<84;t++){
        if(t==0||cmpkey(&ks[t],&ks[t-1])){ if(ks[t].a) circles++; else lines++; if(t){ if(run>maxrun)maxrun=run;} run=1;} else run++;
    }
    if(run>maxrun)maxrun=run;
    if(maxrun==84) return -1; /* all 9 on one block */
    if(print){ printf("SET circles=%d lines=%d pts:",circles,lines); for(int i=0;i<9;i++) printf(" (%d,%d)",px[idx[i]],py[idx[i]]); printf("\n"); }
    return circles;
}
static int idx[9]; static ll total=0;
static void rec(int pos,int start){
    if(pos==9){ int c=count9(idx,0); total++; if(c<0) return; hist[c]++; if(c<=bestc){ if(c<bestc){bestc=c; printf("new best %d\n",c);} count9(idx,1);} return; }
    for(int v=start; v<NP-(8-pos); v++){ idx[pos]=v; rec(pos+1,v+1); }
}
int main(int argc,char**argv){
    W=atoi(argv[1]); H=atoi(argv[2]); int mode=atoi(argv[3]); NP=0;
    for(int x=0;x<W;x++)for(int y=0;y<H;y++){px[NP]=x;py[NP]=y;NP++;}
    bestc = argc>6 ? atoi(argv[6]) : 25; /* print threshold: sets with count <= bestc */
    if(mode==0){ rec(0,0); }
    else { ll S=atoll(argv[4]); srand(atoi(argv[5]));
        for(ll s=0;s<S;s++){ int used[400]={0}; for(int i=0;i<9;i++){int v; do{v=rand()%NP;}while(used[v]); used[v]=1; idx[i]=v;}
            int c=count9(idx,0); total++; if(c<0)continue; hist[c]++; if(c<=bestc){ if(c<bestc){bestc=c; printf("new best %d\n",c);} count9(idx,1);} } }
    printf("total sets %lld; histogram:",total); for(int c=0;c<100;c++) if(hist[c]) printf(" %d:%d",c,hist[c]); printf("\n");
    return 0;
}
