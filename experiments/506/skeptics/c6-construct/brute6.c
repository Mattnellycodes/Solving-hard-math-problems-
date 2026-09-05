/* Exhaustive search over all 6-subsets of a given integer point set for sets with <= MAXC circles.
 * Quadratic form x^2 + c y^2 (c = 1 square lattice, c = 3 triangular lattice in (2a+b, b) coords).
 * Circles identified by exact primitive integer equation vectors (128-bit arithmetic).
 * Input (stdin): N c MAXC, then N lines "x y".  Output: histogram of circle counts and every
 * non-degenerate 6-subset with <= MAXC circles.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
typedef long long ll;
static int N, C, MAXC;
static ll X[512], Y[512], Q[512];
static int *tid;  /* triple id table, index i*N*N + j*N + k, -1 = collinear */
#define TID(i,j,k) tid[((i)*N+(j))*N+(k)]

static i128 gcd128(i128 a, i128 b){ if(a<0)a=-a; if(b<0)b=-b; while(b){i128 t=a%b;a=b;b=t;} return a; }
static i128 det3(i128 a1,i128 a2,i128 a3,i128 b1,i128 b2,i128 b3,i128 c1,i128 c2,i128 c3){
  return a1*(b2*c3-b3*c2)-a2*(b1*c3-b3*c1)+a3*(b1*c2-b2*c1);
}
typedef struct { i128 k[4]; int i,j,l; } key_t_;
static int cmpkey(const void*A,const void*B){
  const key_t_*a=A,*b=B; for(int t=0;t<4;t++){ if(a->k[t]<b->k[t])return -1; if(a->k[t]>b->k[t])return 1;} return 0;
}
int main(){
  if(scanf("%d %d %d",&N,&C,&MAXC)!=3) return 1;
  for(int i=0;i<N;i++){ if(scanf("%lld %lld",&X[i],&Y[i])!=2) return 1; Q[i]=X[i]*X[i]+(ll)C*Y[i]*Y[i]; }
  tid=malloc(sizeof(int)*N*N*N);
  long ntrip=(long)N*(N-1)*(N-2)/6;
  key_t_ *keys=malloc(sizeof(key_t_)*ntrip); long nk=0;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)for(int k=j+1;k<N;k++){
    i128 det=det3(X[i],Y[i],1,X[j],Y[j],1,X[k],Y[k],1);
    if(det==0){ TID(i,j,k)=-1; continue; }
    i128 dn=det3(-Q[i],Y[i],1,-Q[j],Y[j],1,-Q[k],Y[k],1);
    i128 en=det3(X[i],-Q[i],1,X[j],-Q[j],1,X[k],-Q[k],1);
    i128 fn=det3(X[i],Y[i],-Q[i],X[j],Y[j],-Q[j],X[k],Y[k],-Q[k]);
    i128 g=gcd128(gcd128(det,dn),gcd128(en,fn)); det/=g;dn/=g;en/=g;fn/=g;
    if(det<0){det=-det;dn=-dn;en=-en;fn=-fn;}
    keys[nk].k[0]=det;keys[nk].k[1]=dn;keys[nk].k[2]=en;keys[nk].k[3]=fn;keys[nk].i=i;keys[nk].j=j;keys[nk].l=k;nk++;
  }
  qsort(keys,nk,sizeof(key_t_),cmpkey);
  int id=-1;
  for(long t=0;t<nk;t++){ if(t==0||cmpkey(&keys[t],&keys[t-1])!=0) id++; TID(keys[t].i,keys[t].j,keys[t].l)=id; }
  fprintf(stderr,"N=%d triples=%ld distinct circles=%d\n",N,nk,id+1);
  long hist[32]; memset(hist,0,sizeof hist);
  long found=0, total=0;
  int s[6];
  for(s[0]=0;s[0]<N;s[0]++)for(s[1]=s[0]+1;s[1]<N;s[1]++)for(s[2]=s[1]+1;s[2]<N;s[2]++)
  for(s[3]=s[2]+1;s[3]<N;s[3]++)for(s[4]=s[3]+1;s[4]<N;s[4]++)for(s[5]=s[4]+1;s[5]<N;s[5]++){
    int ids[20],m=0,ncol=0;
    for(int a=0;a<6;a++)for(int b=a+1;b<6;b++)for(int c=b+1;c<6;c++){
      int v=TID(s[a],s[b],s[c]); if(v<0){ncol++;continue;}
      int dup=0; for(int t=0;t<m;t++) if(ids[t]==v){dup=1;break;}
      if(!dup) ids[m++]=v;
    }
    total++;
    /* degenerate: all collinear (m==0) or all concyclic (m==1 and ncol==0) */
    if(m==0 || (m==1 && ncol==0)) { hist[31]++; continue; }
    hist[m<31?m:30]++;
    if(m<=MAXC){ found++; printf("circles=%d collinear_triples=%d :",m,ncol); for(int a=0;a<6;a++) printf(" (%lld,%lld)",X[s[a]],Y[s[a]]); printf("\n"); fflush(stdout);}
  }
  fprintf(stderr,"subsets=%ld degenerate=%ld found(<=%d)=%ld\n",total,hist[31],MAXC,found);
  for(int m=1;m<31;m++) if(hist[m]) fprintf(stderr,"  circles=%d : %ld\n",m,hist[m]);
  return 0;
}
