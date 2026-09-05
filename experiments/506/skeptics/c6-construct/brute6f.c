/* Floating-point exhaustive search over 6-subsets of a pool of points for sets with <= MAXC circles.
 * Per triple: circumcentre (cx,cy) and radius r, or collinear flag.  Per 6-subset: cluster the
 * (<=20) circles with a relative tolerance (pairwise), count clusters.  Candidates are printed for
 * exact/independent re-verification.  Input: N MAXC TOL, then N lines "x y" (doubles).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static int N, MAXC; static double TOL;
static double X[400], Y[400];
typedef struct { float cx, cy, r; unsigned char col; } circ;
static circ *tab;
#define T(i,j,k) tab[((i)*N+(j))*N+(k)]
int main(){
  if(scanf("%d %d %lf",&N,&MAXC,&TOL)!=3) return 1;
  for(int i=0;i<N;i++) if(scanf("%lf %lf",&X[i],&Y[i])!=2) return 1;
  tab=malloc(sizeof(circ)*N*N*N);
  double scale=0; for(int i=0;i<N;i++){ if(fabs(X[i])>scale)scale=fabs(X[i]); if(fabs(Y[i])>scale)scale=fabs(Y[i]); }
  long ntrip=0, ncol=0;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)for(int k=j+1;k<N;k++){
    double ax=X[i],ay=Y[i],bx=X[j],by=Y[j],cx=X[k],cy=Y[k];
    double d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
    double lab=hypot(bx-ax,by-ay), lac=hypot(cx-ax,cy-ay), lbc=hypot(cx-bx,cy-by);
    double m=lab; if(lac>m)m=lac; if(lbc>m)m=lbc;
    ntrip++;
    if(fabs(d) < 1e-9*m*m*2){ T(i,j,k).col=1; ncol++; continue; }
    double a2=ax*ax+ay*ay,b2=bx*bx+by*by,c2=cx*cx+cy*cy;
    double ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d, uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d;
    T(i,j,k).cx=ux; T(i,j,k).cy=uy; T(i,j,k).r=hypot(ax-ux,ay-uy); T(i,j,k).col=0;
  }
  fprintf(stderr,"N=%d triples=%ld collinear=%ld scale=%g\n",N,ntrip,ncol,scale);
  long hist[32]; memset(hist,0,sizeof hist); long total=0,found=0;
  int s[6];
  for(s[0]=0;s[0]<N;s[0]++)for(s[1]=s[0]+1;s[1]<N;s[1]++)for(s[2]=s[1]+1;s[2]<N;s[2]++)
  for(s[3]=s[2]+1;s[3]<N;s[3]++)for(s[4]=s[3]+1;s[4]<N;s[4]++)for(s[5]=s[4]+1;s[5]<N;s[5]++){
    circ cs[20]; int m=0, nc=0;
    for(int a=0;a<6;a++)for(int b=a+1;b<6;b++)for(int c=b+1;c<6;c++){
      circ t=T(s[a],s[b],s[c]); if(t.col){nc++;continue;}
      int dup=0;
      for(int q=0;q<m;q++){ double tol=TOL*(1+cs[q].r); if(fabs(cs[q].cx-t.cx)<tol&&fabs(cs[q].cy-t.cy)<tol&&fabs(cs[q].r-t.r)<tol){dup=1;break;} }
      if(!dup) cs[m++]=t;
    }
    total++;
    if(m==0||(m==1&&nc==0)){hist[31]++;continue;}
    hist[m<31?m:30]++;
    if(m<=MAXC){found++; printf("circles=%d collinear=%d :",m,nc); for(int a=0;a<6;a++) printf(" (%.12g,%.12g)",X[s[a]],Y[s[a]]); printf("\n"); fflush(stdout);}
  }
  fprintf(stderr,"subsets=%ld degenerate=%ld found(<=%d)=%ld\n",total,hist[31],MAXC,found);
  for(int m=1;m<31;m++) if(hist[m]) fprintf(stderr,"  circles=%d : %ld\n",m,hist[m]);
  return 0;
}
