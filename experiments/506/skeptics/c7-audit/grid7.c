/* Skeptic c7-audit: exhaustive search over all 7-subsets of a k x k integer grid for the minimum
   number of circles (through >= 3 points) of a non-degenerate (not all collinear / concyclic) set.
   Independent exact arithmetic: a circle through 3 integer points is identified by the primitive
   integer vector (A, B, C, Dn) with  Dn*(x^2+y^2) + A x + B y + C = 0, Dn > 0. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
typedef struct { ll a,b,c,d; } key;
static int cmpkey(const void*x,const void*y){ const key*p=x,*q=y;
  if(p->d!=q->d) return p->d<q->d?-1:1; if(p->a!=q->a) return p->a<q->a?-1:1;
  if(p->b!=q->b) return p->b<q->b?-1:1; if(p->c!=q->c) return p->c<q->c?-1:1; return 0; }
/* returns number of circles; sets *degenerate if all 7 points on one line or one circle */
static int circles7(const int*X,const int*Y,int*degenerate){
  key ks[35]; int nk=0; int collinear=0; key lk[35]; int nl=0;
  for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)for(int k=j+1;k<7;k++){
    ll x1=X[i],y1=Y[i],x2=X[j],y2=Y[j],x3=X[k],y3=Y[k];
    ll A1=x2-x1,B1=y2-y1,R1=-(x2*x2+y2*y2-x1*x1-y1*y1);
    ll A2=x3-x1,B2=y3-y1,R2=-(x3*x3+y3*y3-x1*x1-y1*y1);
    ll det=A1*B2-A2*B1;
    if(det==0){ collinear++;
      /* line key: primitive (A,B,C) with A x + B y = C, sign-normalised */
      ll A=y2-y1,B=x1-x2,C=A*x1+B*y1; ll g=gcdll(gcdll(A,B),C); if(g==0)g=1; A/=g;B/=g;C/=g;
      if(A<0||(A==0&&B<0)){A=-A;B=-B;C=-C;}
      lk[nl].a=A;lk[nl].b=B;lk[nl].c=C;lk[nl].d=0;nl++; continue; }
    ll a=R1*B2-R2*B1, b=A1*R2-A2*R1;           /* a/det, b/det */
    ll c=-(x1*x1+y1*y1)*det - a*x1 - b*y1;     /* c/det */
    ll g=gcdll(gcdll(a,b),gcdll(c,det)); if(det<0)g=-g; a/=g;b/=g;c/=g; ll d=det/g;
    ks[nk].a=a;ks[nk].b=b;ks[nk].c=c;ks[nk].d=d;nk++;
  }
  qsort(ks,nk,sizeof(key),cmpkey);
  int nc=0; for(int i=0;i<nk;i++) if(i==0||cmpkey(&ks[i],&ks[i-1])!=0) nc++;
  *degenerate = (collinear==35) || (nk==35 && nc==1);
  return nc;
}
int main(int argc,char**argv){
  int k=atoi(argv[1]); int N=k*k; int px[64],py[64];
  for(int i=0;i<N;i++){px[i]=i%k;py[i]=i/k;}
  int best=1000; long long hist[40]; memset(hist,0,sizeof hist); int bx[7],by[7]; long long total=0;
  int c[7];
  for(c[0]=0;c[0]<N;c[0]++)for(c[1]=c[0]+1;c[1]<N;c[1]++)for(c[2]=c[1]+1;c[2]<N;c[2]++)
  for(c[3]=c[2]+1;c[3]<N;c[3]++)for(c[4]=c[3]+1;c[4]<N;c[4]++)for(c[5]=c[4]+1;c[5]<N;c[5]++)
  for(c[6]=c[5]+1;c[6]<N;c[6]++){
    int X[7],Y[7]; for(int t=0;t<7;t++){X[t]=px[c[t]];Y[t]=py[c[t]];}
    int deg; int nc=circles7(X,Y,&deg); total++;
    if(deg) continue;
    if(nc<40) hist[nc]++;
    if(nc<best){best=nc; for(int t=0;t<7;t++){bx[t]=X[t];by[t]=Y[t];}
      printf("new min %d: ",nc); for(int t=0;t<7;t++)printf("(%d,%d) ",X[t],Y[t]); printf("\n"); fflush(stdout);}
  }
  printf("grid %dx%d: %lld subsets, min circles (non-degenerate) = %d\n",k,k,total,best);
  for(int i=0;i<40;i++) if(hist[i]) printf("  circles=%d : %lld sets\n",i,hist[i]);
  return 0;
}
