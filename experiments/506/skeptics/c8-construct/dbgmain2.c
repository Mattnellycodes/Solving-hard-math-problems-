#include "grid8_lib.c"
int main(){
  W=5;H=5;N=25; for(int i=0;i<N;i++){px[i]=i%W;py[i]=i/W;}
  int pts[8][2]={{0,2},{1,2},{2,2},{0,3},{2,3},{0,4},{1,4},{2,4}};
  int S[8]; for(int i=0;i<8;i++) S[i]=pts[i][1]*W+pts[i][0];
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
      int isp=0; for(int t=0;t<8;t++) if(fabs(xs[u]-px[S[t]])<1e-9&&fabs(ys[u]-py[S[t]])<1e-9) isp=1; if(isp) continue;
      int deg=0; for(int t=0;t<nb;t++) if(on_block(&B[t],xs[u],ys[u])) deg++;
      if(deg>=3) printf("pt (%g,%g) from blocks %d,%d deg %d\n",xs[u],ys[u],i,j,deg);
      if(deg>best) best=deg;
    }
  }
  printf("best %d\n",best);
}
