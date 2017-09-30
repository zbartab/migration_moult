/*
  migra1: optimal annual routine model for a migratory bird

  developers: Z. Barta       University of Bristol
              R.K. Welham    University of Bristol
	      J.M. McNamara  University of Bristol
	      A.I. Houston   University of Bristol

  maintainer: Z. Barta <zbarta@delfin.unideb.hu>

*/



#include"migra1.h"


extern type_par P[ML];

extern double ******states;

extern FILE *OutPut3;

double G[MT+1][ME+1][ML][N_PATCH];
double oG[MT+1][ME+1][ML][N_PATCH];
double FP[MF+1][ML];
double FOOD[MT+1][ML];
double p_birtok[MT+1][ML];


/* these arrays below introduced because of speed optimization */
double pred_dfe[MF+1][ML];
double fora_dfe[MF+1][ML];
double mifo_dfe[MF+1][ML];
double mipr_dfe[MF+1][ML];
double theta[ME+1][ML];

double p_i[ML][N_PATCH];

double newres(double *, int, int, int, int, int, double, int);
double newfq1(double, int, int);
void iniFP(void);
void inigain(double);
void ReadIn(char *);
void Migration(double *, double *, double *, int, int, int);
double P_birtok(double, double, int);
void iniBinom(void);
void freeP_I(void);
void cp_arr(double *dest, double *src);
void free_arrays(double ******rv);
double ******init_arrays(void);
void free_h_rv(double *****rv);
double *****init_h_rv(void);
void hiba(char *leiras);
double food_dens(double n, int l);


/************************* newfq1 *************************/
double newfq1(double e_exp, int f1, int l)
{
  double nfq1;

  /*  if (MF == 0) return 0.0;*/
  if (f1 == MOL-1) return ((double)MF);	/* finishing moult */
  if (f1 < MOL-1) {    /* growing feathers */
    nfq1 = (((double) f1)+P[l].F.NU/*-P[l].F.foraM*e_exp*/);
    nfq1 = (nfq1 < 0.0 ? 0.0 : nfq1);
    return nfq1;
  }
  if (P[l].F.fixed) nfq1 = ((double) f1) - P[l].F.fora;
  else {
    nfq1 = ((double) f1) - P[l].F.fora*e_exp;
  }
  return (nfq1 < ((double)MOL) ? ((double)MOL) : nfq1);
 }
/************************* end of newfq1 ******************/

/************************* newres *************************/
double newres(double *Pred, int t, int r, int f1, int e, int l, double u, 
	      int action)
{

  double c;
  double absr, a2, u2;

  absr = ((double) r)/TOPR;
  a2 = absr*absr;
  c = P[l].C.base[t] + P[l].C.mass * a2;
  u2 = u*u;
  if ((*Pred = P[l].P.base+u2*(0.1*a2+1.0)*pred_dfe[f1][l]) > 1.0) 
    *Pred = 1.0;
  c += u2*(0.1*a2+1.0)*fora_dfe[f1][l]; 
  if (f1 < MOL) c += P[l].C.moult1;
  if (action==START) c += P[l].C.start;
  else if (action==BIRTOK) c += P[l].C.birtok;
  return c;
}
/************************* end of newres ******************/

#define KOZEL 0.95

/************************* P_birtok *************************/
double P_birtok(double bird, double occ, int l) 
{

  return P[l].L.K*(1.0/(P[l].L.al+(bird/(P[l].L.n_birtok-occ))));
}
/************************* end of P_birtok *************************/


/************************* Migration *************************/
void Migration(double *newr, double *newf1, double *death, int r, int f1, int l)
/* migration costs do depend on reserves and feather quality */
{
  double f1m;
  double absr, a2, expend;

  absr = ((double) r)/TOPR;
  a2 = 0.1*absr*absr+1.0;

  expend = a2*mifo_dfe[f1][l];
  *newr = ((double) r) - expend;
  if(f1>=MOL) {
    f1m = ((double) f1) - expend*P[l].M.feather;
    *newf1 = ( f1m < MOL ? MOL : f1m);
  }
  else {
    f1m = ((double)f1) /*- expend*P[l].M.featherM*/;
    *newf1 = (f1m < 0 ? 0 : f1m);
  }
  *death = a2*mipr_dfe[f1][l];
}
/************************* end of Migration *************************/


/************************* iniBinom *************************/
void iniBinom(void)
{
  int i, l;
  int n = N_PATCH-1;
  int f[N_PATCH];
  double s_p = 0.0;

  f[0] = 1;
  for(i=1; i<N_PATCH; i++) f[i] = i*f[i-1];  
  
  for(l=0; l<ML; l++){
    s_p = 0.0;
    for(i=0; i<N_PATCH; i++) {
      s_p += (p_i[l][i] = 
	      (f[n]/(f[i]*f[n-i]))*pow(P[l].L.p,i)*pow(1.0-P[l].L.p,n-i));
    }
    for(i=0; i<N_PATCH; i++) p_i[l][i] /= s_p;
  }
  
}
/************************* end of iniBinom *************************/



/************************* iniFP *************************/
void iniFP(void)
{
  int i, f1, l;
  double absf, dfe;

  for(l = 0; l < ML; l++) {
    for(i=0; i <=MF; i++){
      if (i < MOL) {
	absf = ((double)(i))/((double)MOL);
      } else {
	absf = ((double) (i-MOL))/TOPF;
      }
      FP[i][l] = P[l].F.worn+pow(absf,P[l].F.flight)*(1.0-P[l].F.worn);
    }
  
    for(f1 = 0; f1 <= MF; f1++) {
      dfe = FP[f1][l];
      /*      dfe = P[l].F.def_coeff*dfe1;*/
      pred_dfe[f1][l] = P[l].P.mass/dfe;
      fora_dfe[f1][l] = P[l].C.propNF*P[l].C.fora + 
	(P[l].C.propF*P[l].C.def)/dfe;
      mifo_dfe[f1][l] = P[l].M.energy/dfe;
      mipr_dfe[f1][l] = P[l].M.preda/dfe;
    }
  }

}
/************************* end of iniFP ******************/

/************************* inigain *************************/
void inigain(double lambda)
{
  double n_fora = 1.0, g;
  int t, e, l, i;

  for(t=0; t<=MT; t++)
    for(l = 0; l < ML; l++){
      for(e=0; e<=ME; e++){
	n_fora = P[l].L.N0;
	g = FOOD[t][l] * theta[e][l];
	for(i = 0; i < N_PATCH; i++) {
	  oG[t][e][l][i] = G[t][e][l][i] =  TOPR * 
	    (g - P[l].L.p*(N_PATCH-1)*P[l].L.omega + i*P[l].L.omega);
	}
      }
      p_birtok[t][l] = P_birtok(1.0,0,l);
    }
}
/************************* end of inigain ******************/

/************************* food_dens *************************/
double food_dens(double n, int l)
{
  return exp(-P[l].L.DD * ((n/P[l].L.N0)-1.0));
}
/************************* food_dens *************************/


/************************************* ReadIn ******************************/
void ReadIn(char *filename)
{
  FILE *InPut;
  char filen[30];
  char *ext = ".ini";
  char junk[50];
  int errors = 0;
  int t, r, fe, a, e, fo, l, lo, n_p;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,ext);
  if ((InPut = fopen(filen,"rt")) == NULL) {
    printf("ERROR: Can't open the input file: %s\n", filen);
    exit (1);
  }
  /* feather dynamics */
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%d", &(P[l].F.fixed))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].F.fora))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].F.flight))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].F.worn))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].F.NU))) errors++;
  /* migration */
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].M.feather))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].M.energy))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].M.preda))) errors++;
  /* metabolisms */
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.fora))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.propF))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.def))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.birtok))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.start))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.keep))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.moult1))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.mass))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.baseM))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].C.baseS))) errors++;
  /* predation */
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].P.base))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].P.mass))) errors++;
  /* food */
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.THETA))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.MAXEPS))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.FOOD))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.p))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.omega))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.n))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.N0))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.DD))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.pe))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.n_birtok))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.K))) errors++;
  if(1!=fscanf(InPut,"%s",junk)) errors++;
  for(l = 0; l < ML; l++) 
    if (1 != fscanf(InPut, "%lf", &(P[l].L.al))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P[0].tol))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P[0].alpha))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P[0].K))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P[0].simit))) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &r)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &t)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &fe)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &fo)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &a)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &e)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &lo)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &n_p)) errors++;
  if (errors > 0) {
    printf("ERROR: in reading input parameters from file %s at %d\n", filen, 
	   errors);
    exit(1);
  }
  if ((MR!=r) || (MF!=fe) || (MT!=t) || (ME!=e) || (MA!=a) || (MOL!=fo) || 
      (ML!=lo) || (N_PATCH!=n_p)) {
    printf("ERROR: the compiled grid not the same as the assumed one!\n");
    printf("\tCompiled:\tMT=%d\tMR=%d\tMF=%d\tME=%d\tMA=%d\tMOL=%d\tML=%d\tNP=%d\n",
	   MT, MR, MF, ME, MA, MOL, ML, N_PATCH);
    printf("\tAssumed:\tMT=%d\tMR=%d\tMF=%d\tME=%d\tMA=%d\tMOL=%d\tML=%d\tNP=%d\n",
	   t, r, fe, e, a, fo, lo, n_p);
    exit(1);
  }
  P[0].n_alpha = 1.0-3.0*P[0].alpha;
  P[0].alpha2 = 2.0*P[0].alpha;
  P[0].K = 1.0/P[0].K;
  for(l = 0; l < ML; l++) {
    P[l].C.fora *= TOPR;
    P[l].C.def *= TOPR;
    P[l].C.start *= TOPR;
    P[l].C.keep *= TOPR;
    P[l].C.birtok *= TOPR;
    P[l].C.moult1 *= TOPR;
    P[l].C.mass *= TOPR;
    P[l].M.energy *= TOPR;
    P[l].F.fora /= TOPR;
    P[l].M.feather /= TOPR;
    P[l].C.propNF = 1.0 - P[l].C.propF;
    P[l].F.NU = 1.0-P[l].F.NU;
    P[l].F.foraM = P[l].F.fora * MOL;
    P[l].M.featherM = P[l].M.feather * MOL;
    P[l].M.feather *= TOPF;
    P[l].F.fora *= TOPF;
    for(e=0; e<=ME; e++) 
      theta[e][l] = pow(((double) P[l].L.THETA), ((double) ME-e));
    for(t=0; t<=MT; t++){
      FOOD[t][l] = P[l].L.FOOD + 
	P[l].L.MAXEPS*sin(2.0*PI*((double) t)/((double) MT) - PI/2.0);
      if(FOOD[t][l] < 0.0) FOOD[t][l] = 0.0;
      P[l].C.base[t] = TOPR*
	(P[l].C.baseM+P[l].C.baseS*sin(2.0*PI*((double) t)/((double) MT) - 
				       PI/2.0));
    }
  }
  fclose(InPut);
}
/************************************* ReadIn ***********************/


/************************* print_action *************************/
void cp_arr(double *dest, double *src)
{
  int i;

  for(i = 0; i < MDEC; i++) dest[i] = src[i];
}
/************************* print_action *************************/

/************************* Dec_Err *************************/
double Dec_Err(double x)
{
  return 1.0/(1.0 + P[0].K*x);
}
/************************* Dec_Err *************************/

/************************* init_arrays *************************/
double ******init_arrays(void)
{
  int t, r, f1, e, a;
  double ******rv;

  rv = (double ******) calloc(MT+1,sizeof(double *****));
  if(rv==NULL) hiba("init_array 1");
  for(t = 0; t <= MT; t++) {
    rv[t] = (double *****) calloc(MR+1,sizeof(double ****));
    if(rv[t]==NULL) hiba("init_array 2");
    for(r = 0; r <=MR; r++) {
      rv[t][r] = (double ****) calloc(MF+1,sizeof(double ***));
      if(rv[t][r] == NULL) hiba("init_array 3");
      for(f1 = 0; f1 <= MF; f1++) {
	rv[t][r][f1] = (double ***) calloc(ME+1,sizeof(double **));
	if(rv[t][r][f1] == NULL) hiba("init_array 4");
	for(e = 0; e <= ME; e++) {
	  rv[t][r][f1][e] = (double **) calloc(MA+1,sizeof(double *));
	  if(rv[t][r][f1][e] == NULL) hiba("init_array 5");
	  for(a = 0; a <= MA; a++) {
	    rv[t][r][f1][e][a] = (double *) calloc(ML,sizeof(double));
	    if(rv[t][r][f1][e][a] == NULL) hiba("init_array 6");
	  }
	}	
      }
    }
  }
  return rv;
}
/************************* init_arrays *************************/

/************************* free_arrays *************************/
void free_arrays(double ******rv)
{
  int t, r, f1, e, a;

  for(t = 0; t <= MT; t++) {
    for(r = 0; r <=MR; r++) {
      for(f1 = 0; f1 <= MF; f1++) {
	for(e = 0; e <= ME; e++) {
	  for(a = 0; a <= MA; a++) {
	    free(rv[t][r][f1][e][a]);
	  }
	  free(rv[t][r][f1][e]);
	}	
	free(rv[t][r][f1]);
      }
      free(rv[t][r]);
    }
    free(rv[t]);
  }
  free(rv);
}
/************************* init_arrays *************************/

/************************* init_h_rv *************************/
double *****init_h_rv(void)
{
  int r, f1, e, a;
  double *****rv;

  rv = (double *****) calloc(MR+1,sizeof(double ****));
  if(rv==NULL) hiba("init_array 2");
  for(r = 0; r <=MR; r++) {
    rv[r] = (double ****) calloc(MF+1,sizeof(double ***));
    if(rv[r] == NULL) hiba("init_array 3");
    for(f1 = 0; f1 <= MF; f1++) {
      rv[r][f1] = (double ***) calloc(ME+1,sizeof(double **));
      if(rv[r][f1] == NULL) hiba("init_array 4");
      for(e = 0; e <= ME; e++) {
	rv[r][f1][e] = (double **) calloc(MA+1,sizeof(double *));
	if(rv[r][f1][e] == NULL) hiba("init_array 5");
	for(a = 0; a <= MA; a++) {
	  rv[r][f1][e][a] = (double *) calloc(ML,sizeof(double));
	  if(rv[r][f1][e][a] == NULL) hiba("init_array 6");
	}
      }	
    }
  }
  
  return rv;
}
/************************* init_h_rv *************************/

/************************* free_h_rv *************************/
void free_h_rv(double *****rv)
{
  int r, f1, e, a;

  for(r = 0; r <=MR; r++) {
    for(f1 = 0; f1 <= MF; f1++) {
      for(e = 0; e <= ME; e++) {
	for(a = 0; a <= MA; a++) {
	  free(rv[r][f1][e][a]);
	}
	free(rv[r][f1][e]);
      }	
      free(rv[r][f1]);
    }
    free(rv[r]);
  }
  free(rv);
}
/************************* free_h_rv *************************/

/************************* hiba *************************/
void hiba(char *leiras)
{
  fprintf(stderr, "ERROR: %s\n", leiras);
  exit(1);
}
/************************* hiba *************************/

