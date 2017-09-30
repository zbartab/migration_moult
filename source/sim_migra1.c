/*
  migra1: optimal annual routine model for a migratory bird

  developers: Z. Barta       University of Bristol
              R.K. Welham    University of Bristol
	      J.M. McNamara  University of Bristol
	      A.I. Houston   University of Bristol

  maintainer: Z. Barta <zbarta@delfin.unideb.hu>

*/


#include<string.h>
#include"migra1.h"

#define NBIRDS 10000
#define NYEARS 5
#define FNBIRD ((float) NBIRDS)

extern double ******youngs;
extern double uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];
extern double G[MT+1][ME+1][ML][N_PATCH];
extern double p_birtok[MT+1][ML];
extern double p_i[ML][N_PATCH];

extern type_par P[ML];

extern FILE *OutPut3;

t_ind Bird[NBIRDS];

extern float randf(float);
extern double newres(double *, int, int, int, int, int, double, int);
extern double newfq1(double, int, int);
extern void Migration(double *, double *, double *, int, int, int);
extern void Treatment(int, t_ind *);
extern void Food_Treatment(int, int);

void ynormalise(double *, int );
void Ind(char *, int);
void MoultStr(int, char *);
void ActionStr(int migact, int act, char *S);
void WriteOut(int, int, int, char *, char *, char *, t_ind, double, FILE *);
int smooth_r(double);
int smooth_f1(double);
int my_rdev(double *, int);

/************************* Ind *************************/
void Ind(char *filename, int tt)
{
  char filen[30], S1[30], S3[30];
  int t, r, f1, e, a, act, migact, ido, l;
  double lambda, u, pmort, e_exp;
  int year, i;
  FILE *OutPut1;
  double rf, ff1, p0, p1, keep;
  int alive;
  int WRITE;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".dat");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  fprintf(OutPut1, "year\tweek\tbird\tevent\tmoult1\taction\tres\tfea1\texper\tbrood\tloc\tuval\n");

  ynormalise(&lambda, SUM);

  for(year = 0; year < NYEARS; year++){
    alive = 0;
    WRITE = year >= 0;
    /*    for(l = 0; l < ML; l++)*/
    for(t = 0; t < MT; t++){
      ido = year*52 + t;
      if (tt != -1) Food_Treatment(ido, t);
      for(i = 0; i < NBIRDS; i++){
	if (year == 0 && Bird[i].birth == t) Bird[i].birth = -1;
	if (Bird[i].alive && Bird[i].birth == -1){
	  migact = 0;
	  keep = 0.0;
	  pmort = 0.0;
	  if (t==0) alive++;
	  if (tt != -1) Treatment(ido, &(Bird[i]));
	  if (Bird[i].r <= 0.0) {
	    if (WRITE) WriteOut(year,t,i,"S",S1,S3,Bird[i],u,OutPut1);
	    Bird[i].alive = 0;
	  }
	  else if (Bird[i].pred) {
	    if (WRITE) WriteOut(year,t,i,"P",S1,S3,Bird[i],u,OutPut1);
	    Bird[i].alive = 0;
	  }
	  else {
	    r = Bird[i].r;
	    f1 = Bird[i].f1;
	    e = Bird[i].e;
	    a = Bird[i].a;
	    l = Bird[i].l;
	    strcpy(S1,"");
	    if (randf(1.0) < moult1[t][r][f1][e][a][l]) {
	      f1 = 0;
	      strcpy(S1, "S");
	    } 
	    if(randf(1.0) < action[t][r][f1][e][a][l][D_FLYNORTH]) {
	      Migration(&rf,&ff1,&pmort,r,f1,l);
	      r = smooth_r(rf);
	      f1 = smooth_f1(ff1);
	      if (l>0) l--;
	      a = 0; 
	      migact = FLYNORTH;
	    } else if(randf(1.0) < action[t][r][f1][e][a][l][D_FLYSOUTH]) {
	      Migration(&rf,&ff1,&pmort,r,f1,l);
	      r = smooth_r(rf);
	      f1 = smooth_f1(ff1);
	      if (l<(ML-1)) l++;
	      a = 0; 
	      migact = FLYSOUTH;
	    } 
	    MoultStr(f1, S1);
	    if (randf(1.0) < pmort) Bird[i].pred = 1;
	    p0 = action[t][r][f1][e][a][l][0];
	    p1 = action[t][r][f1][e][a][l][1];
	    p0 = p0/(p0+p1);
	    if(a == 0) {
	      if(randf(1.0) < p0) {
		act = DELAY;
		u = uvalue[t][r][f1][e][a][l][0];
	      } else {
		act = BIRTOK;
		u = uvalue[t][r][f1][e][a][l][1];
	      }
	    } else if(a == 1) {
	      if(randf(1.0) < p0) {
		act = DELAY;
		u = uvalue[t][r][f1][e][a][l][0];
	      } else {
		act = START;
		u = uvalue[t][r][f1][e][a][l][1];
	      }
	    } else if(a==MA) {
	      if(randf(1.0) < p0) {
		act = ABORT;
		u = uvalue[t][r][f1][e][a][l][0];
	      } else {
		act = ABANDON;
		keep = P[l].C.keep;
		u = uvalue[t][r][f1][e][a][l][1];
	      }
	    } else {
	      if(randf(1.0) < p0) {
		act = ABORT;
		u = uvalue[t][r][f1][e][a][l][0];
	      } else {
		act = KEEP;
		keep = P[l].C.keep;
		u = uvalue[t][r][f1][e][a][l][1];
	      }
	    }
	    ActionStr(migact,act,S3);
	    if (WRITE) WriteOut(year,t,i,"A",S1,S3,Bird[i],u,OutPut1);

	    e_exp = newres(&pmort,t,r,f1,e,l,u,act);
	    ff1 = newfq1(e_exp,f1,l);
	    rf = ((double) r) + u*G[t][e][l][my_rdev(p_i[l],N_PATCH)] - 
				   e_exp - keep; 
	    Bird[i].r = smooth_r(rf);
	    Bird[i].f1 = smooth_f1(ff1);
	    if (randf(1.0) < pmort) Bird[i].pred = 1;
	    if (randf(1.0)<P[l].L.pe) Bird[i].e = (e+1 < ME ? e+1 : ME);
	    Bird[i].a = a;
	    if (act == START) Bird[i].a = 2;
	    else if (act == KEEP)  Bird[i].a = (a+1 > MA ? MA : a+1);
	    else if (act == BIRTOK && randf(1.0) < p_birtok[t][l]) 
	      Bird[i].a = 1;
	    else if (act == ABANDON || act == ABORT) Bird[i].a = 1;
	    Bird[i].l = l;
	  }
	}
      }
    }
    fprintf(OutPut3,"%d, %d\n", year, alive);
  }
  fclose(OutPut1);
}
/************************* end of Ind ******************/

#define SF 4

/************************* rdev *************************/
int my_rdev(double *p, int n_cat) 
{
  double coin, sum = 0.0;
  int i;

  coin = randf(1.0);
  for(i = 0; i < n_cat; i++) {
    sum += p[i];
    if (coin <= sum) return i;
  }
  return (n_cat-1);
}
/************************* end of rdev ******************/

/************************* smooth_r *************************/
/* mimic the smoothing of r */

int smooth_r(double rf) 
{
  double p[SF], pr;
  int r[SF];

  rf = (rf < 0.0 ? 0.0 : rf);
  if (rf >= TOPR) {
    r[1] = r[2] = r[3] = MR;
    r[0] = MR - 1;
    p[2] = p[3] = 0.0;
    p[1] = 1.0 - P[0].alpha;
    p[0] = P[0].alpha;
  }
  else {
    r[1] = (int) rf;
    r[2] = r[1] + 1;
    r[3] = ((r[3]=r[2]+1) > MR ? MR : r[3]);
    r[0] = ((r[0]=r[1]-1) > 0 ? r[0] : 0);
    pr = rf - ((double) r[1]);
    p[0] = P[0].alpha*(1.0-pr);
    p[1] = 1.0 - pr*P[0].n_alpha - P[0].alpha2;
    p[2] = pr*P[0].n_alpha + P[0].alpha;
    p[3] = pr*P[0].alpha;
  }

  return r[my_rdev(p,SF)];
}
/************************* end of smooth_r ******************/

/************************* smooth_f1 *************************/
int smooth_f1(double f1) 
{
  double pf1;
  int f1lo, f1hi;

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
  if (randf(1.0)<pf1) return f1hi;
  else return f1lo;
}
/************************* end of smooth_f1 ******************/

/************************* ynormalise *************************/
void ynormalise(double *lambda, int how)
{
  int t, r, f1, e, a, l;
  double tmp, L1, L2;
  int i;

  *lambda = 0.0;
  for(t = 0; t <= MT; t++)
  for(r = 1; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++)
  for(l = 0; l < ML; l++){
    tmp = youngs[t][r][f1][e][a][l];
    if ((how == MAX) && (*lambda < tmp)) *lambda = tmp;
    if (how == SUM) *lambda += tmp;
  }
  L1 = L2 = 0.0;
  i = 0;
  if (*lambda > 0.0){
    for(t = 0; t <= MT; t++)
    for(r = 1; r <= MR; r++)
    for(f1 = 0; f1 <= MF; f1++)
    for(e = 0; e <= ME; e++)
    for(a = 0; a <= MA; a++)
    for(l = 0; l < ML; l++){
      L2 += (youngs[t][r][f1][e][a][l] /= *lambda);
      while(((((float) i)/FNBIRD) < L2) && (i<NBIRDS)){
	Bird[i].alive = 1;
	Bird[i].pred = 0;
	Bird[i].r = r;
	Bird[i].f1 = f1;
	Bird[i].e = e;
	Bird[i].a = a;
	Bird[i].l = l;
	Bird[i].birth= t;
	i++;
      }
    }
  }
  else {
    fprintf(OutPut3,"ERROR: Something is wrong with lambda in ynormalise\n");
    exit(1);
  }
}
/************************* end of ynormalise ******************/


/************************* MoultStr *************************/
void MoultStr(int f, char *S)
{
  
  if(strlen(S) > 0) return;
  if(f<MOL-1) strcpy(S, "M");
  else if(f==MOL-1) strcpy(S, "E");
  else strcpy(S, "N");
  return;
}
/************************* end of MoultStr ******************/

/************************* MoultStr *************************/
void ActionStr(int migact, int act, char *S)
{
  switch (migact) {
  case FLYSOUTH: strcpy(S, "So"); break;
  case FLYNORTH: strcpy(S, "No"); break;
  default: strcpy(S,""); break;
  }
  switch (act) {
  case DELAY: strcat(S, "D"); break;
  case START: strcat(S, "S"); break;
  case KEEP: strcat(S, "K"); break;
  case ABORT: strcat(S, "AT"); break;
  case ABANDON: strcat(S, "A"); break;
  case BIRTOK: strcat(S, "B"); break;
  default: sprintf(S,"%d", act); break;
  }
}
/************************* end of MoultStr ******************/


/************************* WriteOut *************************/
void WriteOut(int year, int t, int i, char *S0, char *S1, char *S3,
	      t_ind Bird, double u, FILE *OutPut)
     /* this function writes out the bird's state and its decision at
	the begining of week t */
{
  fprintf(OutPut, "%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\n",
	  year, t, i, S0, S1, S3, Bird.r, Bird.f1,
	  Bird.e,Bird.a,Bird.l,u);

}
  /************************* WriteOut *************************/

