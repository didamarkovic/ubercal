#include <stdlib.h>
// #include <fstream>
// #include <iomanip>
#include <stdio.h>
#include <math.h>
// #include <string>

#define NRANSI
#define ITMAX 1000
#define TOL 2.0e-4
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

int ncom;
double *pcom,*xicom;
double (*nrfunc)(double []);
void err_handler(const char*);

const double R=0.61803399;
const double C=(1.0-0.61803399);

// find minimum of f(x)
double golden(double ax, double bx, double cx, double (*f)(double), 
	      double tol, double *xmin)
{
  double f1,f2,x0,x1,x2,x3;
  
  x0=ax;
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax)) { x1=bx; x2=bx+C*(cx-bx); } 
  else                           { x2=bx; x1=bx-C*(bx-ax); }
  f1=(*f)(x1);
  f2=(*f)(x2);
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) { x0=x1; x1=x2; x2=R*x1+C*x3; f1=f2; f2=(*f)(x2); }
    else         { x3=x2; x2=x1; x1=R*x2+C*x0; f2=f1; f1=(*f)(x1); }
  }
  if (f1 < f2) { *xmin=x1; return f1; } 
  else         { *xmin=x2; return f2; }
}

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	    double (*func)(double []))
{
  double *dvector(long,long);
  void   free_dvector(double*,long,long);

  void linmin(double p[], double xi[], int n, double *fret,
	      double (*func)(double []));
  int i,ibig,j;
  double *pt,*ptt,*xit;
  double del,fp,fptt,t;
  
  pt=dvector(1,n);
  ptt=dvector(1,n);
  xit=dvector(1,n);
  *fret=(*func)(p);
  for (j=1;j<=n;j++) pt[j]=p[j];
  for (*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);
      if (fabs(fptt-(*fret)) > del) {
	del=fabs(fptt-(*fret));
	ibig=i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      free_dvector(xit,1,n);
      free_dvector(ptt,1,n);
      free_dvector(pt,1,n);
      return;
    }
    if (*iter == ITMAX) { fprintf(stderr,"powell exceeding maximum iterations.\n"); break; }
    for (j=1;j<=n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
	linmin(p,xit,n,fret,func);
	for (j=1;j<=n;j++) {
	  xi[j][ibig]=xi[j][n];
	  xi[j][n]=xit[j];
	}
      }
    }
  }
}

void linmin(double p[], double xi[], int n, 
	    double *fret, double (*func)(double []))
{
  double *dvector(long,long);
  void   free_dvector(double*,long,long);
  double brent(double ax, double bx, double cx,
	      double (*f)(double), double tol, double *xmin);
  double f1dim(double x);
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	      double *fc, double (*func)(double));
  int j;
  double xx,xmin,bx,ax;
  double fx,fb,fa;

  ncom=n;
  pcom=dvector(1,n);
  xicom=dvector(1,n);
  nrfunc=func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom,1,n);
  free_dvector(pcom,1,n);
}

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	    double *xmin)
{
  int iter;
  double a,b,etemp,tol1,tol2,u,v,w,x,xm;
  double d=0.0,e=0.0;
  double fu,fv,fw,fx,p,q,r;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }
  err_handler("Too many iterations in brent");
  *xmin=x;
  return fx;
}

double f1dim(double x)
{
  double *dvector(long,long);
  void   free_dvector(double*,long,long);

  int j;
  double f;
  double *xt;
  
  xt=dvector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_dvector(xt,1,ncom);
  return f;
}

void mnbrak(double *ax, double *bx, double *cx, 
	    double *fa, double *fb, double *fc, double (*func)(double))
{
  double ulim,u,r,q;
  double fu,dum;
  
  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,(*func)(u))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}
#undef NRANSI
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef TOL
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

