// Exhaustive search: all 7-subsets of a W x H integer grid, count circles through >=3 points exactly.
// Circle key: normalised center (ux/d, uy/d) with gcd removed, plus R = r^2 * d^2 (integer).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a;}
int W,H,N; int px[256],py[256];
typedef struct {ll ux,uy,d,R;} key;
static int count_circles(int *s, int n, int *lines3){
  key ks[64]; int nk=0; *lines3=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++){
    ll x1=px[s[i]],y1=py[s[i]],x2=px[s[j]],y2=py[s[j]],x3=px[s[k]],y3=py[s[k]];
    ll d=2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    if(d==0){(*lines3)++;continue;}
    ll s1=x1*x1+y1*y1,s2=x2*x2+y2*y2,s3=x3*x3+y3*y3;
    ll ux=((s2-s1)*(y3-y1)-(s3-s1)*(y2-y1));
    ll uy=((x2-x1)*(s3-s1)-(x3-x1)*(s2-s1));
    ll g=gcdll(gcdll(ux,uy),d); if(d<0)g=-g; ux/=g;uy/=g;d/=g;
    ll R=(x1*d-ux)*(x1*d-ux)+(y1*d-uy)*(y1*d-uy);
    int found=0; for(int q=0;q<nk;q++) if(ks[q].ux==ux&&ks[q].uy==uy&&ks[q].d==d&&ks[q].R==R){found=1;break;}
    if(!found){ks[nk].ux=ux;ks[nk].uy=uy;ks[nk].d=d;ks[nk].R=R;nk++;}
  }
  return nk;
}
static int degenerate(int *s,int n){ // all on one line or one circle
  int l; int c=count_circles(s,n,&l);
  if(c==0) return 1;
  if(c>1) return 0;
  // c==1: all non-collinear triples on the same circle; check every point on it: if any collinear triple exists
  // with c==1, then some point off the circle would create another circle unless ... just check all points on the circle
  // find the circle from first non-collinear triple
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++){
    ll x1=px[s[i]],y1=py[s[i]],x2=px[s[j]],y2=py[s[j]],x3=px[s[k]],y3=py[s[k]];
    ll d=2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)); if(d==0)continue;
    ll s1=x1*x1+y1*y1,s2=x2*x2+y2*y2,s3=x3*x3+y3*y3;
    ll ux=((s2-s1)*(y3-y1)-(s3-s1)*(y2-y1)); ll uy=((x2-x1)*(s3-s1)-(x3-x1)*(s2-s1));
    ll R=(x1*d-ux)*(x1*d-ux)+(y1*d-uy)*(y1*d-uy);
    for(int q=0;q<n;q++){ ll xx=px[s[q]],yy=py[s[q]]; if((xx*d-ux)*(xx*d-ux)+(yy*d-uy)*(yy*d-uy)!=R) return 0;}
    return 1;
  }
  return 1;
}
int main(int argc,char**argv){
  W=atoi(argv[1]);H=atoi(argv[2]); int K=7; int target=argc>3?atoi(argv[3]):10;
  N=0; for(int y=0;y<H;y++)for(int x=0;x<W;x++){px[N]=x;py[N]=y;N++;}
  int s[8]; int best=1000; long long cnt=0, hits=0;
  // iterate combinations; symmetry: require point s[0] ... no symmetry reduction, brute force
  for(s[0]=0;s[0]<N;s[0]++)for(s[1]=s[0]+1;s[1]<N;s[1]++)for(s[2]=s[1]+1;s[2]<N;s[2]++)
  for(s[3]=s[2]+1;s[3]<N;s[3]++)for(s[4]=s[3]+1;s[4]<N;s[4]++)for(s[5]=s[4]+1;s[5]<N;s[5]++)
  for(s[6]=s[5]+1;s[6]<N;s[6]++){
    cnt++; int l; int c=count_circles(s,K,&l);
    if(c<=target){ if(degenerate(s,K)) continue; hits++;
      if(c<best){best=c; printf("NEW BEST %d circles, %d collinear triples:",c,l); for(int i=0;i<K;i++)printf(" (%d,%d)",px[s[i]],py[s[i]]); printf("\n"); fflush(stdout);}
    }
  }
  printf("grid %dx%d: subsets=%lld, nondegenerate subsets with <=%d circles: %lld, best=%d\n",W,H,cnt,target,hits,best);
  return 0;
}
