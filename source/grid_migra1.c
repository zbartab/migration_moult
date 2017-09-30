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

c-----------------------------------------------------------------------------
c ********************** GRID ROUTINES - MOULT VERSION ***********************
c-----------------------------------------------------------------------------
c  AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c  STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED
c
c-----------------------------------------------------------------------------

rewritten in C by Z. Barta
*/

#include"migra1.h"

extern double ******repval;
extern double ******states;
extern double ******youngs;
extern double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];

extern FILE *OutPut3;

extern type_par P[ML];

void initialise(int, double, int);
void savevalues_back(void);
void savevalues_for(double);
void setup_n_for(double);
double gridread(double *****, double, double, int, int, int);
double gridget(int, int, int, int, int, int, int);
void normalise(double *, double *, int, int);
void gridwrite(double *****, double, double, int, int, int, double);
void gridinc(double *****, int, int, int, int, int, double);
void disturb_for(double lambda);

/************************* gridinc *************************/
void gridinc(double *****strv, int r, int f1, int e, int a, int l,
	     double inc)
{

  if (r == 0) return;
  strv[r][f1][e][a][l] += inc;
}
/************************* end of gridinc ******************/


/************************* gridwrite *************************/
void gridwrite(double *****strv, double r, double f1, int e, 
	       int a, int l, double inc)
{
  double pr, pf1, qr, qf1;
  int f1lo, f1hi;
  double p1, p2, p3, p4;
  int r1, r2, r3, r4;

  if (inc <= 0.0 || r <= 0.0) return;
  else if (r >= TOPR) {
    r2 = r3 = r4 = MR;
    r1 = MR - 1;
    p3 = p4 = 0.0;
    p2 = 1.0 - P[0].alpha;
    p1 = P[0].alpha;
  }
  else {
    r2 = (int) r;
    r3 = r2 + 1;
    r4 = ((r4=r3+1) > MR ? MR : r4);
    r1 = ((r1=r2-1) > 0 ? r1 : 0);
    pr = r - ((double) r2);
    p1 = P[0].alpha*(1.0-pr);
    p2 = 1.0 - pr*P[0].n_alpha - P[0].alpha2;
    p3 = pr*P[0].n_alpha + P[0].alpha;
    p4 = pr*P[0].alpha;
  }
  if (f1 <= (MOL-0.9)) {
    f1lo = (int) f1;
    f1hi = f1lo+1;
    pf1 = f1-((double) f1lo);
  }
  else if (f1 <= MOL) {
    f1lo = MOL;
    f1hi = MOL;
    pf1 = 1.0;
  }
  else if (f1 >= MF) {
    f1lo = MF;
    f1hi = MF;
    pf1 = 0.0;
  }
  else {
    f1lo = (int) f1;
    f1hi = f1lo + 1;
    pf1 = f1 - ((double) f1lo);
  }
  qr = 1.0 - pr;
  qf1 = 1.0 - pf1;

  if(r2!=0) {
    gridinc(strv,r1,f1lo,e,a,l,p1*qf1*inc);
    gridinc(strv,r1,f1hi,e,a,l,p1*pf1*inc);
    gridinc(strv,r2,f1lo,e,a,l,p2*qf1*inc);
    gridinc(strv,r2,f1hi,e,a,l,p2*pf1*inc);
  }

  gridinc(strv,r3,f1lo,e,a,l,p3*qf1*inc);
  gridinc(strv,r3,f1hi,e,a,l,p3*pf1*inc);
  gridinc(strv,r4,f1lo,e,a,l,p4*qf1*inc);
  gridinc(strv,r4,f1hi,e,a,l,p4*pf1*inc);
}
/************************* end of gridwrite ******************/


/************************* normalise *************************/
void normalise(double *lambda, double *diff, int how, int act)
{
  int r, f1, e, a, l;
  double tmp, tmp2;

  *lambda = 0.0;
  *diff = 0.0;
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    if (how == MAX) {
      tmp = repval[0][r][f1][e][a][l];
      if (*lambda < tmp) *lambda = tmp;
    }
    if (how == SUM) {
      *lambda += states[0][r][f1][e][a][l];
    }
    }
  /*  if(how == MAX)   *lambda = repval[0][MR][MF][ME][1][0];*/
  if (*lambda > 0.0){
    tmp = 0.0;
    for(r = 0; r <= MR; r++)
    for(f1 = 0; f1 <= MF; f1++)
    for(l = 0; l < ML; l++)
    for(e = 0; e <= ME; e++)
    for(a = 0; a <= MA; a++){
      if (how == MAX) {
	if (act==NORM) repval[0][r][f1][e][a][l] /= *lambda;
	*diff += fabs(repval[0][r][f1][e][a][l] -
		      repval[MT][r][f1][e][a][l]);
      } else {
	if (act==NORM) states[0][r][f1][e][a][l] /= *lambda;
	*diff += (tmp2 = fabs(states[0][r][f1][e][a][l] -
		      states[MT][r][f1][e][a][l]));
	if (tmp < tmp2) tmp = tmp2;
      }
    }
  }
  else {
    fprintf(OutPut3,
	    "ERROR: Something is wrong with lambda in normalise: %f\n", 
	    *lambda);
  }
}
/************************* end of normalise ******************/


/************************* gridread *************************/
double gridread(double *****rv, double r, double f1, int e, int a, int l)
{
  double pr, pf1, qr, qf1;
  double gr = 0.0;
  double p1, p2, p3, p4;
  int r1, r2, r3, r4;
  int f1lo, f1hi;

  if (r <= 0.0) return 0.0;
  else if (r >= TOPR) {
    r2 = r3 = r4 = MR;
    r1 = MR - 1;
    p3 = p4 = 0.0;
    p2 = 1.0 - P[0].alpha;
    p1 = P[0].alpha;
  }
  else {
    r2 = (int) r;
    r3 = r2 + 1;
    r4 = ((r4=r3+1) > MR ? MR : r4);
    r1 = ((r1=r2-1) > 0 ? r1 : 0);
    pr = r - ((double) r2);
    p1 = P[0].alpha*(1.0-pr);
    p2 = 1.0 - pr*P[0].n_alpha - P[0].alpha2;
    p3 = pr*P[0].n_alpha + P[0].alpha;
    p4 = pr*P[0].alpha;
  }
  if (f1 <= MOL-0.9) {
    f1lo = (int) f1;
    f1hi = f1lo+1;
    pf1 = f1-((double) f1lo);
  }
  else if (f1 <= MOL) {
    f1lo = MOL;
    f1hi = MOL;
    pf1 = 1.0;
  }
  else if (f1 >= ((double)MF)) {
    f1lo = MF;
    f1hi = MF;
    pf1 = 0.0;
  }
  else {
   f1lo = (int) f1;
   f1hi = f1lo + 1;
   pf1 = f1 - ((double) f1lo);
  }
  qr = 1.0 - pr;
  qf1 = 1.0 - pf1;

  if(r2!=0) {
    gr = qf1*(p1*rv[r1][f1lo][e][a][l]  + 
	      p2*rv[r2][f1lo][e][a][l]) + 
         pf1*(p1*rv[r1][f1hi][e][a][l]  + 
	      p2*rv[r2][f1hi][e][a][l]);
  }

  gr += qf1*(p3*rv[r3][f1lo][e][a][l]  + 
	     p4*rv[r4][f1lo][e][a][l]) + 
        pf1*(p3*rv[r3][f1hi][e][a][l]  + 
	     p4*rv[r4][f1hi][e][a][l]);
  return gr;
}
/************************* end of gridread ******************/

/************************* savevalues_back *************************/
void savevalues_back(void)
{
  int r, f1, e, a, l;

  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    repval[MT][r][f1][e][a][l] = 
      repval[0][r][f1][e][a][l];
  }
}
/************************* end of savevalues_back ******************/


/************************* savevalues_for *************************/
void savevalues_for(double lambda)
{
  int r, f1, e, a, l;

  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    states[MT][r][f1][e][a][l] = 
      states[0][r][f1][e][a][l];
  }
}
/************************* end of savevalues_for ******************/

/************************* setup_n_for *************************/
void setup_n_for(double lambda)
{
  int r, f1, e, a, l;
  double sum_n;

  sum_n = 0.0;
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    sum_n += states[0][r][f1][e][a][l];
  }
  lambda /= sum_n;
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    states[0][r][f1][e][a][l] *= lambda;
    states[MT][r][f1][e][a][l] = 
      states[0][r][f1][e][a][l];
  }
}
/************************* setup_n_for ******************/



/************************* initialise *************************/
void initialise(int t, double initval, int what)
{
  int r, f1, e, a, l;

  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(a = 0; a <= MA; a++){
    if (what == BACKWARD) {
      if (r == 0) repval[t][r][f1][e][a][l] = 0.0;
      else repval[t][r][f1][e][a][l] = initval;
      action[MT-1][r][f1][e][a][l][D_DELAY] = 1.0;
    } else {
      states[t][r][f1][e][a][l] = 0.0;
    }
  }
  for(l=0; l<ML; l++) states[t][MR][MF][ME][0][l] = initval;
}
/************************* end of initialise ******************/


