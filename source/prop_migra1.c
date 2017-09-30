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

c------------------------------------------------------------------------------
c ***** BACKWARDS AND FORWARDS VALUE PROPAGATION - MOULT VERSION **************
c------------------------------------------------------------------------------
c    AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c    STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED
c
c------------------------------------------------------------------------------

rewritten in C by Z. Barta
*/

#include"migra1.h"



extern double ******repval;
extern double ******states;
extern double ******youngs;
extern double *****h_rv;


extern double uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];
extern double theta[ME+1][ML];
extern double FOOD[MT+1][ML];
extern double G[MT+1][ME+1][ML][N_PATCH];
extern double p_birtok[MT+1][ML];

extern double p_i[ML][N_PATCH];


extern type_par P[ML];

extern void initialise(int, double, int);
extern double newres(double *, int, int, int, int, int, double, int);
extern double newfq1(double, int, int);
extern double gridread(double *****, double, double, int, int, int);
extern void gridwrite(double *****, double, double, int, int, int, double);
extern double P_birtok(double, double, int);
extern void Migration(double *, double *, double *, int, int, int);
extern void copy_states(int t);
extern void dumpm_states(double *****st);
extern double food_dens(double n, int l);

double backvalue(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l,
		 double u, int actn, double *pb1, double keep);
void stepforwards(int, int, int, double);
void sf_i(int how, int t, int t1, int r, int f1, int e, int e1, int a, 
	  int a1, int l, int actn, double keep, double u, 
	  double *pb1, double w);
double sum_states(int t);


/************************* stepforwards ***********************/
void stepforwards(int t, int how, int code, double lambda)
{
  double w, pmort, rm, f1m;
  int r, f1, e, e1, a, l, t1, i;
  double N_for[ME+1], n_fora, N_birtok;
  double N_bird;
  double p_fs, p_fn;
  double p_delay, g;
  double pb1[4];

#if DEBUG
  printf("1: %d %f\n",t, sum_states(t));
  dumpm_states(states[t]);
  dumpm_states(h_rv);
#endif
  t1 = (t+1) % MT;
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(r = MR; r >= 0; r--)
  for(f1 = MF; f1 >= 0; f1--)
  for(a = 0; a <= MA; a++){
    h_rv[r][f1][e][a][l] = 0.0;
  }
  for(l = 0; l < ML; l++)
  for(e = 0; e <= ME; e++)
  for(r = MR; r >= 0; r--)
  for(f1 = MF; f1 >= 0; f1--)
  for(a = 0; a <= MA; a++){
    if(f1 >= MOL) {
      p_fn=moult1[t][r][f1][e][a][l];
      w = states[t][r][f1][e][a][l];
      states[t][r][0][e][a][l] += w * p_fn;
      states[t][r][f1][e][a][l] -= w * p_fn;
    }
    w = states[t][r][f1][e][a][l]; 
    /*    if (a<2) {*/
		p_fs = p_fn = 0.0;
		if (l < (ML-1)) {
			p_fs = action[t][r][f1][e][a][l][D_FLYSOUTH];
			Migration(&rm,&f1m,&pmort,r,f1,l);
			gridwrite(h_rv,rm,f1m,e,0,l+1,(1.0-pmort)*p_fs*w);
		} 
		if (l > 0) {
			p_fn = action[t][r][f1][e][a][l][D_FLYNORTH];
			Migration(&rm,&f1m,&pmort,r,f1,l);
			gridwrite(h_rv,rm,f1m,e,0,l-1,(1.0-pmort)*p_fn*w);
		} 
		h_rv[r][f1][e][a][l] += states[t][r][f1][e][a][l]*(1.0 - p_fs - p_fn);
      /*    } else {
      h_rv[r][f1][e][a][l] += states[t][r][f1][e][a][l]; 
      }*/
	}

  for(l = 0; l < ML; l++){
    N_birtok = 0.0;
    n_fora = 0.0;
    N_bird = 0.0;
    for(e = 0; e <= ME; e++){
      N_for[e] = 0.0;
      for(r = 0; r <= MR; r++)
	for(f1 = 0; f1 <= MF; f1++)
	  for(a = 0; a <= MA; a++){
	    w = h_rv[r][f1][e][a][l]; 
	    if(w>0.0) {
	      p_delay = action[t][r][f1][e][a][l][D_DELAY];
	      N_for[e] += w*(p_delay*uvalue[t][r][f1][e][a][l][D_DELAY] + 
			     (1.0-p_delay)*uvalue[t][r][f1][e][a][l][D_BIRTOK]);
	      if(a>0) N_birtok += w;
	      else {
		N_bird += w*(1.0-p_delay);
	      }
	    }
	  }
      n_fora += N_for[e] * theta[e][l];
    }
    p_birtok[t][l] = P_birtok(N_bird,N_birtok,l);
    for(e = 0; e <= ME; e++) {
      g = (FOOD[t][l]*food_dens(n_fora,l)) * theta[e][l];
      for(i = 0; i < N_PATCH; i++) {
	G[t][e][l][i] =  TOPR * 
	  (g - P[l].L.p*(N_PATCH-1)*P[l].L.omega + i*P[l].L.omega);
      }
    }
  }
#if DEBUG
  printf("2: %d %f\n", t, sum_states(t));
  dumpm_states(states[t]);
  dumpm_states(h_rv);
#endif
  if (code == INIT) initialise(t1, 0.0, FORWARD);
  for(a = 0; a <= MA; a++)
  for(f1 = 0; f1 <= MF; f1++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++)
  for(l = 0; l < ML; l++){
    w = h_rv[r][f1][e][a][l];
    if (w > 0.0) {
      e1 = (e < ME ? e+1 : ME);
      pb1[0] = (e==ME ? 1.0 : (1.0-P[l].L.pe));
      pb1[1] = (e==ME ? 0.0 : P[l].L.pe);
      pb1[2] = pb1[3] = 0.0;
      
      if(a == 0) {		/* territory search */
	sf_i(how,t,t1,r,f1,e,e1,0,0,l,DELAY,0.0,
	     uvalue[t][r][f1][e][a][l][D_DELAY],
	     pb1,w*action[t][r][f1][e][a][l][D_DELAY]);
	pb1[0] = (e==ME ? p_birtok[t][l] : (1.0-P[l].L.pe)*p_birtok[t][l]);
	pb1[1] = (e==ME ? 0.0 : P[l].L.pe*p_birtok[t][l]);
	pb1[2] = (e==ME ? (1.0-p_birtok[t][l]) : (1.0-P[l].L.pe)*(1.0-p_birtok[t][l]));
	pb1[3] = (e==ME ? 0.0 : P[l].L.pe*(1.0-p_birtok[t][l]));
	sf_i(how,t,t1,r,f1,e,e1,0,1,l,BIRTOK,0.0,
	     uvalue[t][r][f1][e][a][l][D_BIRTOK],
	     pb1,w*action[t][r][f1][e][a][l][D_BIRTOK]);
      } else if (a == 1) {	/* start a brood */
	sf_i(how,t,t1,r,f1,e,e1,1,1,l,DELAY,0.0,
	     uvalue[t][r][f1][e][a][l][D_DELAY],
	     pb1,w*action[t][r][f1][e][a][l][D_DELAY]);
	sf_i(how,t,t1,r,f1,e,e1,2,2,l,START,0.0,
	     uvalue[t][r][f1][e][a][l][D_START],
	     pb1,w*action[t][r][f1][e][a][l][D_START]);
      } else if (a >= MA) {	/* abandon a brood */
	sf_i(how,t,t1,r,f1,e,e1,1,1,l,ABORT,0.0,
	     uvalue[t][r][f1][e][a][l][D_ABORT],
	     pb1,w*action[t][r][f1][e][a][l][D_ABORT]);
	sf_i(how,t,t1,r,f1,e,e1,1,1,l,ABANDON,P[l].C.keep,
	     uvalue[t][r][f1][e][a][l][D_ABANDON],
	     pb1,w*action[t][r][f1][e][a][l][D_ABANDON]);
      } else {		/* keep a brood */
	sf_i(how,t,t1,r,f1,e,e1,1,1,l,ABORT,0.0,
	     uvalue[t][r][f1][e][a][l][D_ABORT],
	     pb1,w*action[t][r][f1][e][a][l][D_ABORT]);
	sf_i(how,t,t1,r,f1,e,e1,a+1,a+1,l,KEEP,P[l].C.keep,
	     uvalue[t][r][f1][e][a][l][D_KEEP],
	     pb1,w*action[t][r][f1][e][a][l][D_KEEP]);
      }
    }
  }
}
/************************* end of stepforwards ******************/

/************************* sf_i *************************/
void sf_i(int how, int t, int t1, int r, int f1, int e, int e1, int a, 
	  int a1, int l, int actn, double keep, double u, 
	  double *pb1, double w) 
{
  double e_exp, newf1, newr, pmort;
  int i;

  e_exp = newres(&pmort,t,r,f1,e,l,u,actn);
  newf1 = newfq1(e_exp,f1,l);
  w *= (1.0-pmort);
  if (u>0.0) {
    for(i=0; i<N_PATCH; i++) {
      newr = ((double) r) + u*G[t][e][l][i] - keep - e_exp;
      gridwrite(states[t1],newr,newf1,e, a1,l,p_i[l][i]*pb1[0]*w);
      if(pb1[1]>0.0) gridwrite(states[t1],newr,newf1,e1,a1,l,p_i[l][i]*pb1[1]*w);
      if(pb1[2]>0.0) gridwrite(states[t1],newr,newf1,e,  a,l,p_i[l][i]*pb1[2]*w);
      if(pb1[3]>0.0) gridwrite(states[t1],newr,newf1,e1, a,l,p_i[l][i]*pb1[3]*w);
    }
  } else {
    newr = ((double) r) - keep - e_exp;
    gridwrite(states[t1],newr,newf1,e, a1,l,pb1[0]*w);
    if(pb1[1]>0.0) gridwrite(states[t1],newr,newf1,e1,a1,l,pb1[1]*w);
    if(pb1[2]>0.0) gridwrite(states[t1],newr,newf1,e,  a,l,pb1[2]*w);
    if(pb1[3]>0.0) gridwrite(states[t1],newr,newf1,e1, a,l,pb1[3]*w);
  }
  if (actn == ABANDON && how == ALL) {
    gridwrite(states[t1],0.5*TOPR,((double)MF),0,0,l,P[l].L.n*w);
    gridwrite(youngs[t1],0.5*TOPR,((double)MF),0,0,l,P[l].L.n*w);
  }
}
/************************* sf_i *************************/

/************************* backvalue *************************/
double backvalue(int t, int t1, int r, int f1, int e, int e1, int a, int a1, 
		 int l, double u, int actn, double *pb1, double keep)
{
  double newr, newf1, pmort;
  double bv = 0.0;
  double e_exp;
  int i;

  e_exp = newres(&pmort,t,r,f1,e,l,u,actn);
  newf1 = newfq1(e_exp,f1,l);
  if (u > 0.0) {
    for(i = 0; i < N_PATCH; i++) {
      newr = ((double) r) +  u*G[t][e][l][i] - keep - e_exp;
      bv += p_i[l][i]*pb1[0]*gridread(repval[t1], newr, newf1, e,  a1, l);
      if(pb1[1]>0.0)bv+=p_i[l][i]*pb1[1]*gridread(repval[t1],newr,newf1,e1,a1,l);
      if(pb1[2]>0.0)bv +=p_i[l][i]*pb1[2]*gridread(repval[t1],newr,newf1,e,a, l);
      if(pb1[3]>0.0)bv +=p_i[l][i]*pb1[3]*gridread(repval[t1],newr,newf1,e1,a,l);
    }
  } else {
    newr = ((double) r) - keep - e_exp; 
    bv = pb1[0]*gridread(repval[t1], newr, newf1, e,  a1, l);
    if(pb1[1]>0.0) bv += pb1[1]*gridread(repval[t1], newr, newf1, e1,  a1, l);
    if(pb1[2]>0.0) bv += pb1[2]*gridread(repval[t1], newr, newf1, e,   a, l);
    if(pb1[3]>0.0) bv += pb1[3]*gridread(repval[t1], newr, newf1, e1,  a, l);
  }
  if (actn == ABANDON) {
    bv += P[l].L.n*gridread(repval[t1],0.5*TOPR, ((double)MF), 0, 0, l);
  }
  return bv*(1.0-pmort);
}
/************************* end of backvalue ******************/


/************************* sum_states *************************/
double sum_states(int t) 
{
  int r, f1, e, a, l;
  double s = 0.0;

  for(a = 0; a <= MA; a++)
  for(f1 = 0; f1 <= MF; f1++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++)
  for(l = 0; l < ML; l++) {
    s += states[t][r][f1][e][a][l];
  }
  return s;
}
/************************* sum_states *************************/

