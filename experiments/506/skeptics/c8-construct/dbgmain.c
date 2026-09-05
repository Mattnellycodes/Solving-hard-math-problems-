#include "grid8_lib.c"
int main(){
  W=5;H=5;N=25; for(int i=0;i<N;i++){px[i]=i%W;py[i]=i/W;}
  int pts[8][2]={{0,2},{1,2},{2,2},{0,3},{2,3},{0,4},{1,4},{2,4}};
  int S[8]; for(int i=0;i<8;i++) S[i]=pts[i][1]*W+pts[i][0];
  /* replicate moebius_best with prints */
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
    printf("block %d isline %d mask %x cx %g cy %g r %g\n",nb-1,b.isline,b.mask,b.cx,b.cy,b.r);
  }
  int nbb; int r=moebius_best(S,&nbb); printf("moebius_best=%d nb=%d\n",r,nbb);
  return 0;
}
