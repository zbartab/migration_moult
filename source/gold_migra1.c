/*
  migra1: optimal annual routine model for a migratory bird

  developers: Z. Barta       University of Bristol
              R.K. Welham    University of Bristol
	      J.M. McNamara  University of Bristol
	      A.I. Houston   University of Bristol

  maintainer: Z. Barta <zbarta@delfin.unideb.hu>

*/

/*

original version:

c--------------------------------------------------------------------------
C**** GOLDEN RATIO SEARCH FOR BEST U VALUE - MOULT/MIGRATION VERSION ******
c--------------------------------------------------------------------------
c  AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL		OCTOBER 1999
c  STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED
c
c--------------------------------------------------------------------------


rewritten in C by Z. Barta
*/

#include"migra1.h"


extern double backvalue(int t, int t1, int r, int f1, int e, int e1, int a, 
			int a1, int l, double u, int actn, double *pb1, 
			double keep);

extern int Sz;

void maxuv(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	   double *pb1, double keep, int actn, double uleft, double uright, 
	   double *umax, double *vmax);
void golden(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	    double *pb1, double keep, int actn, double ul, double um, double ur, 
	    double vm, double *umax, double *vmax);
void csucs(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	   double *pb1, double keep, int actn, double uleft, double uright, 
	   double *umax, double *vmax, int besti, double bestv, double *ugrid);

#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define GR 0.61803399
#define GC 0.38196601

/************************* golden *************************/
void golden(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	    double *pb1, double keep, int actn, double ul, double um, double ur, 
	    double vm, double *umax, double *vmax)
{
  double u0, u1, u2, u3;
  double v1, v2;

  u0 = ul;
  u3 = ur;
  if ((ur-um) > (um-ul)){
    u1 = um;
    v1 = vm;
    u2 = um+GC*(ur-um);
    v2 = backvalue(t,t1,r,f1,e,e1,a,a1,l,u2,actn,pb1,keep);
  } else {
    u2 = um;
    v2 = vm;
    u1 = um-GC*(um-ul);
    v1 = backvalue(t,t1,r,f1,e,e1,a,a1,l,u1,actn,pb1,keep);
  }
  
  while ((u3-u0) > UTOL /* *(u1+u2)*/ ) {
    if (v2 > v1) {
      SHFT3(u0,u1,u2,GR*u1+GC*u3);
      v1 = v2;
      v2 = backvalue(t,t1,r,f1,e,e1,a,a1,l,u2,actn,pb1,keep);
    } else {
      SHFT3(u3,u2,u1,GR*u2+GC*u0);
      v2 = v1;
      v1 = backvalue(t,t1,r,f1,e,e1,a,a1,l,u1,actn,pb1,keep);
    }
  }
  if (v1>v2) {
    *umax = u1;
    *vmax = v1;
  } else {
    *umax = u2;
    *vmax = v2;
  }
  return;
}
/************************* end of golden ******************/



/************************* maxuv *************************/
void maxuv(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	   double *pb1, double keep, int actn, double uleft, double uright, 
	   double *umax, double *vmax)
{
  int i;
  double du, u, v;
  double *ugrid;

  *umax = 0.0;
  *vmax = 0.0;
  if (r <= 0) return;
  if ((ugrid = (double *) calloc(UGRAIN+1,sizeof(double))) == NULL) {
    printf("ERROR: Not enough memory: maxuv 1\n");
    exit(1);
  }
  du = (uright-uleft)/((double) UGRAIN);
  for(i = 0; i <= UGRAIN; i++){
    u = uleft + ((double) i)*du;
    ugrid[i] = backvalue(t,t1,r,f1,e,e1,a,a1,l,u,actn,pb1,keep);
  }
  if(ugrid[0]>=ugrid[1]) {
    csucs(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,uleft,uright,&u,&v,0,ugrid[0],ugrid);
    if(v > *vmax) {*vmax = v; *umax = u;}
  }
  for(i = 1; i < UGRAIN; i++) {	
    if(ugrid[i]>ugrid[i-1] && ugrid[i]>=ugrid[i+1]) {
      csucs(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,uleft,uright,&u,&v,i,ugrid[i],ugrid);
      if(v > *vmax) {*vmax = v; *umax = u;}
    }
  }
  if(ugrid[UGRAIN]>ugrid[UGRAIN-1]) {
    csucs(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,uleft,uright,&u,&v,UGRAIN,ugrid[UGRAIN],ugrid);
    if(v > *vmax) {*vmax = v; *umax = u;}
  }
  free(ugrid);
}
/************************* end of maxuv ******************/




/************************* csucs *************************/
void csucs(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, 
	   double *pb1, double keep, int actn, double uleft, double uright, 
	   double *umax, double *vmax, int besti, double bestv, double *ugrid) 
{
  double ul, ur, du, um, vm;
  
  du = (uright-uleft)/((double) UGRAIN);
  if (besti == 0) {
    ul = uleft;
    ur = uleft+du;
    *umax = ul;
    *vmax = ugrid[0];
    do {
      um = 0.5*(ul+ur);
      vm = backvalue(t,t1,r,f1,e,e1,a,a1,l,um,actn,pb1,keep);
      if (vm > ugrid[0]) {
	golden(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,ul,um,ur,vm,umax,vmax);
	break;
      }
      ur = um;
    } while  ((ur-ul) > UTOL);
  }
  else if (besti == UGRAIN) {
    ur = uright;
    ul = uright - du;
    *umax = ur;
    *vmax = ugrid[UGRAIN];
    do {
      um = 0.5*(ul+ur);
      vm = backvalue(t,t1,r,f1,e,e1,a,a1,l,um,actn,pb1,keep);
      if (vm > ugrid[UGRAIN]) {
	golden(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,ul,um,ur,vm,umax,vmax);
	break;
      }
      ul = um;
    } while ((ur-ul)>UTOL);
  }
  else {
    um = uleft + ((double) besti)*du;
    golden(t,t1,r,f1,e,e1,a,a1,l,pb1,keep,actn,um-du,um,um+du,bestv,umax,vmax);
  }
}
/************************* end of csucs ******************/

