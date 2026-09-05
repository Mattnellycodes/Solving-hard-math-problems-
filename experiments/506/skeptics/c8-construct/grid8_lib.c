/* Exhaustive search over 8-subsets of a WxH integer grid for point sets with few circles.
   Euclidean count = number of distinct circles through >=3 points (collinear triples: none).
   For subsets with few blocks, also compute the Moebius optimum: |B| - max degree of any point
   in the plane (intersection of two blocks), i.e. the best count after inverting about that point.
   Usage: ./grid8 W H [euclid_report_threshold] [moebius_threshold]
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct { long a,b,c,d; } key_t_;
static long gcdl(long a,long b){ if(a<0)a=-a; if(b<0)b=-b; while(b){long t=a%b;a=b;b=t;} return a; }
static int W,H,N;
static int px[100],py[100];
/* triple -> circle id, -1 if collinear */
static int *tripid; /* index i*N*N + j*N + k */
static int ncirc=0;
static key_t_ *keys; static int keycap=0;
static int find_key(key_t_ k){
  for(int i=0;i<ncirc;i++) if(keys[i].a==k.a&&keys[i].b==k.b&&keys[i].c==k.c&&keys[i].d==k.d) return i;
  if(ncirc>=keycap){ keycap=keycap?keycap*2:1024; keys=realloc(keys,keycap*sizeof(key_t_)); }
  keys[ncirc]=k; return ncirc++;
}
static long det3(long m[3][3]){
  return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}
static int circle_id(int i,int j,int k){
  long rows[3][4]; int id[3]={i,j,k};
  for(int t=0;t<3;t++){ long x=px[id[t]],y=py[id[t]]; rows[t][0]=x*x+y*y; rows[t][1]=x; rows[t][2]=y; rows[t][3]=1; }
  long coef[4];
  for(int c=0;c<4;c++){ long m[3][3]; for(int r=0;r<3;r++){int cc=0; for(int q=0;q<4;q++) if(q!=c) m[r][cc++]=rows[r][q];} coef[c]=((c&1)?-1:1)*det3(m); }
  if(coef[0]==0) return -1;
  long g=gcdl(gcdl(coef[0],coef[1]),gcdl(coef[2],coef[3])); if(coef[0]<0) g=-g;
  key_t_ kk={coef[0]/g,coef[1]/g,coef[2]/g,coef[3]/g};
  return find_key(kk);
}
/* Moebius evaluation in floating point */
typedef struct { int isline; double cx,cy,r; double ax,ay,bx,by; unsigned mask; } blk;
static double EPS=1e-7;
static int on_block(blk*b,double x,double y){
  if(b->isline){ double dx=b->bx-b->ax,dy=b->by-b->ay; double v=(x-b->ax)*dy-(y-b->ay)*dx; return fabs(v)<EPS*(1+fabs(dx)+fabs(dy))*(1+fabs(x)+fabs(y)); }
  double d=hypot(x-b->cx,y-b->cy); return fabs(d-b->r)<EPS*(1+b->r);
}
static int moebius_best(int *S, int *nblocks_out){
  /* build blocks from the 8 points */
  blk B[64]; int nb=0; unsigned done[64]; int nd=0;
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++){
    unsigned mk=(1u<<i)|(1u<<j)|(1u<<k);
    int seen=0; for(int t=0;t<nd;t++) if((done[t]&mk)==mk){seen=1;break;} if(seen) continue;
    double x1=px[S[i]],y1=py[S[i]],x2=px[S[j]],y2=py[S[j]],x3=px[S[k]],y3=py[S[k]];
    double d=2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    blk b; b.mask=mk;
    if(d==0){ b.isline=1; b.ax=x1;b.ay=y1;b.bx=x2;b.by=y2; }
    else { b.isline=0; double a2=x1*x1+y1*y1,b2=x2*x2+y2*y2,c2=x3*x3+y3*y3;
      b.cx=(a2*(y2-y3)+b2*(y3-y1)+c2*(y1-y2))/d; b.cy=(a2*(x3-x2)+b2*(x1-x3)+c2*(x2-x1))/d; b.r=hypot(x1-b.cx,y1-b.cy); }
    for(int t=0;t<8;t++) if(!(mk>>t&1)) if(on_block(&b,px[S[t]],py[S[t]])) b.mask|=1u<<t;
    done[nd++]=b.mask; B[nb++]=b;
  }
  *nblocks_out=nb;
  /* candidate points: pairwise intersections */
  int best=0;
  for(int i=0;i<nb;i++)for(int j=i+1;j<nb;j++){
    double xs[2],ys[2]; int np=0;
    blk *p=&B[i],*q=&B[j];
    if(p->isline&&q->isline){
      double d1x=p->bx-p->ax,d1y=p->by-p->ay,d2x=q->bx-q->ax,d2y=q->by-q->ay; double den=d1x*d2y-d1y*d2x; if(fabs(den)<1e-12) continue;
      double t=((q->ax-p->ax)*d2y-(q->ay-p->ay)*d2x)/den; xs[0]=p->ax+t*d1x; ys[0]=p->ay+t*d1y; np=1;
    } else if(p->isline||q->isline){
      blk *L=p->isline?p:q, *C=p->isline?q:p;
      double dx=L->bx-L->ax,dy=L->by-L->ay; double len=hypot(dx,dy); dx/=len;dy/=len;
      double fx=L->ax-C->cx,fy=L->ay-C->cy; double t0=-(fx*dx+fy*dy); double hx=fx+t0*dx,hy=fy+t0*dy; double h2=hx*hx+hy*hy; double disc=C->r*C->r-h2;
      if(disc<-1e-9) continue; if(disc<0) disc=0; double s=sqrt(disc);
      xs[0]=L->ax+(t0+s)*dx; ys[0]=L->ay+(t0+s)*dy; xs[1]=L->ax+(t0-s)*dx; ys[1]=L->ay+(t0-s)*dy; np=2;
    } else {
      double dx=q->cx-p->cx,dy=q->cy-p->cy; double d=hypot(dx,dy); if(d<1e-12) continue;
      double a=(p->r*p->r-q->r*q->r+d*d)/(2*d); double h2=p->r*p->r-a*a; if(h2<-1e-9) continue; if(h2<0)h2=0; double h=sqrt(h2);
      double mx=p->cx+a*dx/d,my=p->cy+a*dy/d; xs[0]=mx+h*dy/d; ys[0]=my-h*dx/d; xs[1]=mx-h*dy/d; ys[1]=my+h*dx/d; np=2;
    }
    for(int u=0;u<np;u++){
      /* skip if it is one of the 8 points */
      int isp=0; for(int t=0;t<8;t++) if(fabs(xs[u]-px[S[t]])<1e-9&&fabs(ys[u]-py[S[t]])<1e-9) isp=1; if(isp) continue;
      int deg=0; for(int t=0;t<nb;t++) if(on_block(&B[t],xs[u],ys[u])) deg++;
      if(deg>best) best=deg;
    }
  }
  return nb-best;
}
int main_orig(int argc,char**argv){
  W=atoi(argv[1]); H=atoi(argv[2]); int thr=argc>3?atoi(argv[3]):17; int mthr=argc>4?atoi(argv[4]):17;
  N=W*H; for(int i=0;i<N;i++){px[i]=i%W;py[i]=i/W;}
  tripid=malloc(sizeof(int)*N*N*N);
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)for(int k=j+1;k<N;k++) tripid[(i*N+j)*N+k]=circle_id(i,j,k);
  fprintf(stderr,"grid %dx%d: %d circles from triples\n",W,H,ncirc);
  int *stamp=calloc(ncirc+1,sizeof(int)); int cur=0;
  int S[8]; long long cnt=0; long hist[60]; memset(hist,0,sizeof hist); long mhist[60]; memset(mhist,0,sizeof mhist);
  int minE=999,minM=999;
  /* maximum external degree is at most 7 (orchard), so Moebius count >= nblocks-7; nblocks <= circles + lines (<=... ) */
  for(S[0]=0;S[0]<N;S[0]++)for(S[1]=S[0]+1;S[1]<N;S[1]++)for(S[2]=S[1]+1;S[2]<N;S[2]++)for(S[3]=S[2]+1;S[3]<N;S[3]++)
  for(S[4]=S[3]+1;S[4]<N;S[4]++)for(S[5]=S[4]+1;S[5]<N;S[5]++)for(S[6]=S[5]+1;S[6]<N;S[6]++)for(S[7]=S[6]+1;S[7]<N;S[7]++){
    cnt++; cur++; int c=0, col=0;
    for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++){
      int id=tripid[(S[i]*N+S[j])*N+S[k]];
      if(id<0){col++;continue;}
      if(stamp[id]!=cur){stamp[id]=cur;c++;}
    }
    if(c==0) continue; /* all collinear */
    /* all concyclic check: c==1 and col==0 */
    if(c==1&&col==0) continue;
    hist[c<59?c:59]++;
    if(c<minE) minE=c;
    if(c<=thr){ printf("EUCLID %d :",c); for(int i=0;i<8;i++) printf(" (%d,%d)",px[S[i]],py[S[i]]); printf("\n"); fflush(stdout);}
    /* Moebius: number of blocks = c + (#lines with >=3 pts). lines <= 8 by packing, ext degree <=7 -> only worthwhile if c+lines-7 <= mthr, i.e. c <= mthr+7 (lines>=0) */
    if(c<=mthr+7){
      int nb; int mb=moebius_best(S,&nb);
      if(mb<minM) minM=mb;
      mhist[mb<59?mb:59]++;
      if(mb<=mthr){ printf("MOEBIUS %d (blocks %d, euclid %d):",mb,nb,c); for(int i=0;i<8;i++) printf(" (%d,%d)",px[S[i]],py[S[i]]); printf("\n"); fflush(stdout);}
    }
  }
  printf("grid %dx%d subsets %lld minEuclid %d minMoebius(among c<=%d) %d\n",W,H,cnt,minE,mthr+7,minM);
  printf("euclid hist:"); for(int i=0;i<60;i++) if(hist[i]) printf(" %d:%ld",i,hist[i]); printf("\n");
  printf("moebius hist:"); for(int i=0;i<60;i++) if(mhist[i]) printf(" %d:%ld",i,mhist[i]); printf("\n");
  return 0;
}
