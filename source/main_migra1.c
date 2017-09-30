/*
  migra1: optimal annual routine model for a migratory bird

  developers: Z. Barta       University of Bristol
              R.K. Welham    University of Bristol
	      J.M. McNamara  University of Bristol
	      A.I. Houston   University of Bristol

  maintainer: Z. Barta <zbarta@delfin.unideb.hu>

*/

/*

Original version: 

c------------------------------------------------------------------------------
c***************** MAIN PROGRAM - MOULT/ MIGRATION VERSION ********************
c------------------------------------------------------------------------------
c
c AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED

rewritten in C by Z. Barta <zbarta@delfin.klte.hu>
*/


#include"migra1.h"


#define MAXN 100
#define MAX_BREED 10

#define CONV_BACK 2000
#define CONV_FOR 500

double ucritcl[MT][ME+1][ML];
double unocare[MR+1][MF+1][ME+1][ML];
double vnocare[MR+1][MF+1][ME+1][ML];

double ******repval;
double ******states;
double ******youngs;

double *****h_rv;		/* temporary storage because of migration */


double uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];

double o_uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double o_action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double o_moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];

extern double G[MT+1][ME+1][ML][N_PATCH];
extern double oG[MT+1][ME+1][ML][N_PATCH];
extern double theta[ME+1][ML];
extern double FOOD[MT+1][ML];
extern double p_birtok[MT+1][ML];

type_par P[ML];

long idum=(-13);

int Breed[MAX_BREED][2][ML];
int n_t = 0;

FILE *OutPut3;

extern void maxuv(int t, int t1, int r, int f1, int e, int e1, int a, int a1, int l, double *pb1, double keep,int actn, 
	   double uleft, double uright, double *umax, double *vmax);
extern void initialise(int, double, int);
extern void savevalues_back(void);
extern void savevalues_for(double);
extern void normalise(double *, double *, int, int);
extern void stepforwards(int, int, int, double);
extern void iniFP(void);
extern void inigain(double);
extern void Ind(char *, int);
extern void ReadIn(char *);
extern void Migration(double *, double *, double *, int, int, int);
extern double gridread(double *****, double, double, int, int, int);
extern void iniBinom(void);
extern double Dec_Err(double x);
extern double check_action(int t);
extern void free_arrays(double ******rv);
extern double ******init_arrays(void);
extern void free_h_rv(double *****rv);
extern double *****init_h_rv(void);
extern void hiba(char *);
void disturb_for(double lambda);

void maxV(int, int, int, int, int, int, double *, double *, 
	  double *);
void Backward(void);
void Forward(double *, double *);
void old_Forward(void);
void OutPut(char *);
void Cohort(char *);
void inityoung(double, int);
void PolicyOut(char *);
void BinWrite(char *);
void Treatment(int, t_ind *); 
void Food_Treatment(int, int);
void get_treat(int, char **, char *, char *, int *);
void init_action(void);
void copy_states(int t);
void dumpm_states(double *****st);
void copy_policy(void);
void damp_policy(void);

#define VISSZA 10

/************************* main *************************/
int main(int argc, char *argv[])
{
  char filename[128];
  char filen[128];
  char ext[128];
  double lambda_back, diff_back, d_lambs;
  double lambda_for = 1.0, diff_for;
  double prev3lam, lambda_back2[VISSZA];
  int I, l, i, v;
  int Kapcsolo = 0;
  double simit;
  double n;
  int Sz_f, Sz_b;

  idum = (long) time(NULL);
  if (idum > 0) idum = -idum;

  repval = init_arrays();
  states = init_arrays();
  youngs = init_arrays();
  h_rv = init_h_rv();

  get_treat(argc, argv, filename, ext, &n_t);
  printf("filename: %s\n", filename);
  strcpy(filen,""); strcat(filen,filename); 
  if(n_t > 0) {strcat(filen,"-"); strcat(filen,ext); }
  strcat(filen,".sum");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file\n");
    exit (1);
  }
  fprintf(OutPut3,"Seed: %ld\n", idum);
  for(l=0; l<ML; l++) for(i=0; i<n_t; i++) 
    fprintf(OutPut3, "-b %d %d %d\n", l, Breed[i][0][l], Breed[i][1][l]);
  ReadIn(filename);
  if(n_t>0) {strcat(filename,"-");strcat(filename,ext);}
  iniFP();
  iniBinom();
  initialise(0, 1.0, BACKWARD);
  initialise(0, 1.0, FORWARD);
  init_action();
  lambda_back = 100.0;
  Sz_b = Sz_f = 0;  
  inigain(1.0);
  simit = P[0].simit;
  P[0].simit = 1.0;
  I = 0;
  /*  K_ori = P[0].K;
      P[0].K = 100.0;*/
  do {
    Sz_b++;
    Backward(); 
    if (Kapcsolo==1) BinWrite(filename);
    normalise(&lambda_back,&diff_back,MAX,NORM);
    printf("backpass done, sz= %d lambda= %f diff = %f n= %f\n", 
	   Sz_b, lambda_back, diff_back, P[0].L.n); 
    if(Sz_b>CONV_BACK){
      fprintf(OutPut3, "ERROR: no convergence in first Backward, Sz= %d\n", 
	      Sz_b);
      break;
    }
  } while (diff_back > P[0].tol);
  copy_policy();
  do {
    Sz_f++;
    I++;
    if ((I % 17) == 0) prev3lam = diff_for;
    Forward(&lambda_back,&diff_for);
    if (Kapcsolo==1) BinWrite(filename);
    normalise(&lambda_for, &diff_for, SUM,NONORM);
    printf("forwardpass done, Sz=%d, lambda=%f, diff=%f, n=%f\n", 
	   Sz_f, lambda_for, diff_for, P[0].L.n);
    if(fabs(prev3lam-diff_for)< 0.1*P[0].tol) {
      fprintf(OutPut3, "WARNING: cycling, Sz= %d\n", Sz_f);
      break;
    }
    if(I>CONV_FOR){
      I = 0;
      fprintf(OutPut3, "ERROR: no convergence in first Forward, Sz= %d\n",Sz_f);
      break;
    }
  } while (diff_for > P[0].tol);
  /*  P[0].simit = simit;*/
  n = -1.0;
  v = 0;
  do {
    do{
      Sz_b++;
      n += 1.0;
      P[0].simit = simit/(1.0+(n/1000.0));
      Backward(); 
      if (Kapcsolo==1) BinWrite(filename);
      normalise(&lambda_back,&diff_back,MAX,NORM);
      printf("backpass done, sz= %d lambda= %f diff = %f simit= %f\n", 
	     Sz_b, lambda_back, diff_back, P[0].simit); 
      if(Sz_b>CONV_BACK){
	fprintf(OutPut3, "ERROR: no convergence in Backward, Sz= %d\n", Sz_b);
	break;
      }
    }while(diff_back > P[0].tol);
    damp_policy();
    I = 0;
    do {
      Sz_f++;
      I++;
      if ((I % 17) == 0) prev3lam = diff_for;
      Forward(&lambda_back,&diff_for);
      if (Kapcsolo==1) BinWrite(filename);
      normalise(&lambda_for, &diff_for, SUM,NONORM);
      printf("forwardpass done, Sz=%d, lambda=%f, diff=%f, n=%f\n", 
	     Sz_f, lambda_for, diff_for, P[0].L.n);
      if(fabs(prev3lam-diff_for)< 0.1*P[0].tol) {
	fprintf(OutPut3, "WARNING: cycling, Sz= %d\n", Sz_f);
	break;
      }
      if(I>CONV_FOR){
	fprintf(OutPut3, "ERROR: no convergence in Forward, Sz= %d\n", Sz_f);
	I = 0;
	break;
      }
    } while (diff_for > P[0].tol);
    memmove((lambda_back2+1),lambda_back2,(VISSZA-1)*sizeof(double));
    lambda_back2[0] = lambda_back;
    v++;
    d_lambs = 0.0;
    for(i = 0; i < (VISSZA-1); i++) {
      d_lambs += fabs(lambda_back2[i+1] - lambda_back2[i]);
    }
	d_lambs /= ((double) VISSZA);
    if(d_lambs < P[0].tol) {
      /*fprintf(OutPut3, "WARNING: K is set to a new value, Sz= %d\n", Sz_b);*/
      if(fabs(1.0-lambda_back) < P[0].tol) break;
      /*      else if(P[0].K < 100.0) break;
      else {
	sprintf(filen,"%s%f",filename,P[0].K);
	BinWrite(filen);
	Sz_b = 0;
	P[0].K *= 0.8;
	}*/
    }
    if((Sz_b > CONV_BACK) || 
       ((lambda_for < P[0].tol/1000.0) && (lambda_back < 1.0))) {
      fprintf(OutPut3, "ERROR: no convergence, Sz= %d\n", Sz_b);
      break;
    }
  } while (1) /*(fabs(1.0-lambda_back) > P[0].tol || diff_back > P[0].tol)*/;
  old_Forward();
  BinWrite(filename);
  Ind(filename, -1);
  free_h_rv(h_rv);
  free_arrays(repval);
  free_arrays(states);
  free_arrays(youngs);
  fclose(OutPut3);
  return 0;
}
/************************* end of main ******************/

/************************* get_treat *************************/
void get_treat(int argc, char *argv[], char *filename, char *ext, int *n) 
     /* args format: -b location start end  filename extension */
{
  int l, i = 1, p = 0;
  int name_read = NO;

  strcpy(filename,"baseline");
  strcpy(ext,"con");
  while (i < argc) {
    if (strcmp(argv[i],"-b") == 0) {
      l = atoi(argv[++i]); /* breeding location */
      Breed[p][0][l] = atoi(argv[++i]); /* breeding start */
      Breed[p][1][l] = atoi(argv[++i]); /* breeding finish */
      p++;
    } else {
      if (name_read == NO) {
	strcpy(filename,argv[i]);
	name_read = YES;
      } else {
	strcpy(ext,argv[i]);
	break;
      }
    }
    i++;
  }
  *n = p;
}
/************************* end of get_treat ******************/


/************************* Backward *************************/
void Backward(void)
{
  int t, e, a, r, f1, l, i;
  double vmax, rm, f1m, pmo, sum_vi;
  double v1, v2, vS, vN, ui, pi;

  for(t = 0; t < MT; t++)
    for(e = 0; e <= ME; e++)
      for(l = 0; l < ML; l++){
	for(i=0; i< N_PATCH; i++){
	oG[t][e][l][i] = G[t][e][l][i] = 
	  oG[t][e][l][i] * (1.0-P[0].simit) + G[t][e][l][i] * P[0].simit;
	}
	ucritcl[t][e][l] = P[l].C.keep/G[t][e][l][(N_PATCH-1)/2];
      }
  savevalues_back();
  for(t = MT-1; t >= 0; t--){
  for(a = 0; a <= MA; a++)
  for(f1 = 0; f1 <= MF; f1++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++)
  for(l = 0; l < ML; l++) {
    maxV(t,r,f1,e,a,l,&v1, uvalue[t][r][f1][e][a][l],
	 action[t][r][f1][e][a][l]);
    h_rv[r][f1][e][a][l] = v1;
  }
  for(l = 0; l < ML; l++) 
  for(e = 0; e <= ME; e++)
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(a = 0; a <= MA; a++){
    /*    if(a < 2) {*/
      vN = vS = 0.0;
      if (l > 0) {		/* migrate north */
	Migration(&rm, &f1m, &pmo, r, f1, l);
	vN = (1.0-pmo)*gridread(h_rv,rm,f1m,e,0,l-1);
      } 
      if (l < (ML-1)) {		/* migrate south */
	Migration(&rm, &f1m, &pmo, r, f1, l);
	vS = (1.0-pmo)*gridread(h_rv,rm,f1m,e,0,l+1);
      }
      v1 = h_rv[r][f1][e][a][l];
      if(v1 > vS && v1 > vN) 
	vmax = v1;
      else if(vS > vN)
	vmax = vS;
      else vmax = vN;
      if(vmax<=0.0) {
	  action[t][r][f1][e][a][l][D_FLYNORTH] = 
	    action[t][r][f1][e][a][l][D_FLYSOUTH] = 
	    repval[t][r][f1][e][a][l] = 0.0;
	  
      } else {
	action[t][r][f1][e][a][l][D_FLYSOUTH] = 
	  Dec_Err((vmax-vS)/vmax);
	action[t][r][f1][e][a][l][D_FLYNORTH] = 
	  Dec_Err((vmax-vN)/vmax);
	rm = Dec_Err((vmax-v1)/vmax);
	if(l==0) action[t][r][f1][e][a][l][D_FLYNORTH] = 0.0;
	if(l==(ML-1)) action[t][r][f1][e][a][l][D_FLYSOUTH] = 0.0;
	sum_vi = action[t][r][f1][e][a][l][D_FLYNORTH] + 
	  action[t][r][f1][e][a][l][D_FLYSOUTH] + rm;
	action[t][r][f1][e][a][l][D_FLYSOUTH] /= sum_vi;
	action[t][r][f1][e][a][l][D_FLYNORTH] /= sum_vi;
	rm /= sum_vi;
	repval[t][r][f1][e][a][l] = v1 * rm + 
	  vN * action[t][r][f1][e][a][l][D_FLYNORTH] +
	  vS * action[t][r][f1][e][a][l][D_FLYSOUTH];
      }
      /*   } else {
      repval[t][r][f1][e][a][l] = h_rv[r][f1][e][a][l];
      }   */ 
    if(f1 >= MOL) {
      v2 = repval[t][r][0][e][a][l];
      v1 = repval[t][r][f1][e][a][l];
      vmax = (v1 > v2 ? v1 : v2);
      if(vmax <= 0.0) {
	pi = moult1[t][r][f1][e][a][l] = 0.0;
      } else {
	ui = Dec_Err((vmax-v2)/vmax);
	pi = moult1[t][r][f1][e][a][l] = ui/(ui + Dec_Err((vmax-v1)/vmax));
      }
      repval[t][r][f1][e][a][l] = pi*v2 + (1.0-pi)*v1; 
    } else {
      moult1[t][r][f1][e][a][l] = 0.0;
    }
  }
  }
}
/************************* end of Backward *****************/


/************************* Forward *************************/
void Forward(double *lambda, double *diff)
{
  int t;

  savevalues_for(*lambda);
  for(t = 0; t < MT; t++) stepforwards(t, ALL, INIT, 1.0);
}
/************************* end of Forward ******************/

/************************* old_Forward *************************/
void old_Forward(void)
{
  int t;
  double lambda, diff;

  inityoung(0.0, INIT);
  savevalues_for(1.0);
  for(t = 0; t < MT; t++) stepforwards(t, ALL, INIT, 1.0);
  normalise(&lambda, &diff, SUM,NONORM);
  fprintf(OutPut3,"forwards convergence, lambda=%f, diff=%f, n=%f\n", 
	  lambda, diff, P[0].L.n);
  
}
/************************* end of old_Forward ******************/

/************************* inityoung *************************/
void inityoung(double initval, int code)
{
  int t, r, f1, e, a, l;
  
  for(t = 0; t <= MT; t++)
  for(r = 0; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(e = 0; e <= ME; e++)
  for(l = 0; l < ML; l++)
  for(a = 0; a <= MA; a++){
    if (code == INIT) {
      if (r == 0) youngs[t][r][f1][e][a][l] = 0.0;
      else youngs[t][r][f1][e][a][l] = initval;
    }
    else {
      states[t][r][f1][e][a][l] = 
	youngs[t][r][f1][e][a][l];
    }
  }
}
/************************* end of inityoung ******************/

/************************* BinWrite *************************/
void BinWrite(char *filename)
{
  char filen[30];
  FILE *BinOut;
  int n_y, n_a, n_m;
  int t, r, f1, e, a, l;

  n_y = (MT+1)*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML;
  n_m = MT*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML;
  n_a = MT*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML*MDEC;
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".food");
  if ((BinOut = fopen(filen,"w")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  for(t = 0; t <= MT; t++) {
    for(e = 0; e <= ME; e++) {
      for(l = 0; l < ML; l++) {
	fprintf(BinOut,"%f\n",G[t][e][l][(N_PATCH-1)/2]);
      }
    }
  }
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".rvb");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fwrite(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  for(t = 0; t <= MT; t++) {
    for(r = 0; r <=MR; r++) {
      for(f1 = 0; f1 <= MF; f1++) {
	for(e = 0; e <= ME; e++) {
	  for(a = 0; a <= MA; a++) {
	    for(l = 0; l < ML; l++) {
	      if (fwrite(&(repval[t][r][f1][e][a][l]), sizeof(double),1,BinOut)
		  != 1)
		hiba("BinWrite 2a");
	    }
	  }
	}
      }
    }
  }
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".stb");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fwrite(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  for(t = 0; t <= MT; t++) {
    /*  t = 0;*/
  for(r = 0; r <=MR; r++) {
    for(f1 = 0; f1 <= MF; f1++) {
      for(e = 0; e <= ME; e++) {
	for(a = 0; a <= MA; a++) {
	  for(l = 0; l < ML; l++) {
	    if (fwrite(&(states[t][r][f1][e][a][l]), sizeof(double),1,BinOut) 
		!= 1) 
	      hiba("BinWrite 2b");
	  }
	}
      }
    }
  }
  }
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".uvb");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fwrite(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinWrite 3");
  if (fwrite(uvalue, sizeof(double),n_a,BinOut) != n_a) hiba("BinWrite 4");
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".acb");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fwrite(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  if (fwrite(action, sizeof(double),n_a,BinOut) != n_a) hiba("BinWrite 5");
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".mob");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fwrite(&n_m,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  if (fwrite(moult1, sizeof(double),n_m,BinOut) != n_m) hiba("BinWrite 6");
  fclose(BinOut);
}
/************************* BinWrite *************************/

/************************* Treatment *************************/
void Treatment(int t, t_ind *Bird)
{
  /* dump function for simulated treatments,
     working version in `moult_treat.c' */
}
/************************* Treatment *************************/

/************************* Food_Treatment *************************/
void Food_Treatment(int ido, int t)
{
  /* dump function for simulated treatments, 
     working version in `moult_treat.c' */
}
/************************* Food_Treatment *************************/

/************************* maxV *************************/
void maxV(int t, int r, int f1, int e, int a, int l,
	  double *vmax, double *umax, double *acmax)
{
  double vi[MDEC];
  double pb1[4];
  double vm, sum_vi;
  int i, t1, e1;
  
  *vmax = 0.0;
  for(i = 0; i < MDEC; i++)  {
    vi[i] = umax[i] = acmax[i] = 0.0;
  }
  if (r <= 0) {
    acmax[D_DELAY] = 1.0;		/* if r==0 default action DELAY */
    return;
  }
  t1 = (t+1) % MT;

  e1 = (e < ME ? e+1 : ME);
  pb1[0] = (e==ME ? 1.0 : (1.0-P[l].L.pe));
  pb1[1] = (e==ME ? 0.0 : P[l].L.pe);
  pb1[2] = pb1[3] = 0.0;
  if(a==0) {			/* get territory */
    maxuv(t,t1,r,f1,e,e1,0,0,l,pb1,0.0,DELAY,0.0,1.0,umax,vi); /* delay */
    pb1[0] = (e==ME ? p_birtok[t][l] : (1.0-P[l].L.pe)*p_birtok[t][l]);
    pb1[1] = (e==ME ? 0.0 : P[l].L.pe*p_birtok[t][l]);
    pb1[2] = (e==ME ? (1.0-p_birtok[t][l]) : (1.0-P[l].L.pe)*(1.0-p_birtok[t][l]));
    pb1[3] = (e==ME ? 0.0 : P[l].L.pe*(1.0-p_birtok[t][l]));
    maxuv(t,t1,r,f1,e,e1,0,1,l,pb1,0.0,BIRTOK,0.0,1.0,(umax+D_BIRTOK),(vi+D_BIRTOK));
  } else if(a==1) {			/* start to breed */
    maxuv(t,t1,r,f1,e,e1,1,1,l,pb1,0.0,DELAY,0.0,1.0,umax,vi); /* delay */
    unocare[r][f1][e][l] = umax[D_DELAY];
    vnocare[r][f1][e][l] = vi[D_DELAY];
    maxuv(t,t1,r,f1,e,e1,2,2,l,pb1,0.0,START,0.0,1.0,(umax+D_START),(vi+D_START));
  } else if(a >= MA) {
    umax[D_ABORT] = unocare[r][f1][e][l];
    vi[D_ABORT] = vnocare[r][f1][e][l];
    if (ucritcl[t][e][l] <= 1.0){
      maxuv(t,t1,r,f1,e,e1,1,1,l,pb1,P[l].C.keep,ABANDON,ucritcl[t][e][l],1.0,(umax+D_ABANDON),
	    (vi+D_ABANDON));
    }
  } else {
    umax[D_ABORT] = unocare[r][f1][e][l];
    vi[D_ABORT] = vnocare[r][f1][e][l];
    if (ucritcl[t][e][l] <= 1.0){
      maxuv(t,t1,r,f1,e,e1,a+1,a+1,l,pb1,P[l].C.keep,KEEP,ucritcl[t][e][l],1.0,(umax+D_KEEP),
	    (vi+D_KEEP));
    }
  }  
  sum_vi = 0.0;
  vm = (vi[D_DELAY] > vi[D_BIRTOK] ? vi[D_DELAY] : vi[D_BIRTOK]);
  if(vm<=0.0) {
    acmax[D_ABORT] = 1.0;
    *vmax = 0.0;
  } else {
    if(a>1 && ucritcl[t][e][l] > 1.0){
      *vmax = vi[D_ABORT];
      acmax[D_ABORT] = 1.0;
      acmax[D_KEEP] = 0.0;
      umax[D_KEEP] = 0.0;
    } else {
    sum_vi = (acmax[D_ABORT] = Dec_Err((vm-vi[D_ABORT])/vm)) +
      (acmax[D_KEEP] = Dec_Err((vm-vi[D_KEEP])/vm));
    *vmax = (acmax[D_ABORT] /= sum_vi) * vi[D_ABORT] +
      (acmax[D_KEEP] /= sum_vi) *vi[D_KEEP];
    }
  }
  return;
}
/************************* maxV *************************/

/************************* init_action *************************/
void init_action(void) 
{
  int t, r, f1, e, a, l;

  for(t = MT-1; t >= 0; t--)
  for(a = 0; a <= MA; a++)
  for(f1 = 0; f1 <= MF; f1++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++)
  for(l = 0; l < ML; l++) {
    action[t][r][f1][e][a][l][D_DELAY] = 1.0;
  }
}
/************************* init_action *************************/

/************************* copy_rv *************************/
void copy_states(int t)
{
  int r, f1, e, a, l;

  for(r = 0; r <=MR; r++) {
    for(f1 = 0; f1 <= MF; f1++) {
      for(e = 0; e <= ME; e++) {
	for(a = 0; a <= MA; a++) {
	  for(l = 0; l < ML; l++) {
	    h_rv[r][f1][e][a][l] = states[t][r][f1][e][a][l];
	  }
	}
      }	
    }
  }
}
/************************* copy_rv *************************/

/************************* copy_rv *************************/
void dumpm_states(double *****st)
{
  int r, f1, e, a, l;
  FILE *dump;

  dump = fopen("dump.txt","w+");
  if(dump==NULL) hiba("dump_states");


  for(r = 0; r <=MR; r++) {
    for(f1 = 0; f1 <= MF; f1++) {
      for(e = 0; e <= ME; e++) {
	for(a = 0; a <= MA; a++) {
	  for(l = 0; l < ML; l++) {
	    fprintf(dump,"%d %d %d %d %d %20.15f\n",r,f1,e,a,l,
		   st[r][f1][e][a][l]);
	  }
	}
      }	
    }
  }
  fclose(dump);
}
/************************* copy_rv *************************/

/************************* copy_policy *************************/
void copy_policy(void)
{
  int t, r, f1, e, a, l, d;

  for(t = 0; t < MT; t++) {
    for(r = 0; r <=MR; r++) {
      for(f1 = 0; f1 <= MF; f1++) {
	for(e = 0; e <= ME; e++) {
	  for(a = 0; a <= MA; a++) {
	    for(l = 0; l < ML; l++) {
	      for(d = 0; d < MDEC; d++) {
		o_uvalue[t][r][f1][e][a][l][d] = uvalue[t][r][f1][e][a][l][d];
		o_action[t][r][f1][e][a][l][d] = action[t][r][f1][e][a][l][d];
	      }
	      o_moult1[t][r][f1][e][a][l] = moult1[t][r][f1][e][a][l];

	    }
	  }
	}
      }
    }
  }

}
/************************* copy_policy *************************/

/************************* damp_policy *************************/
void damp_policy(void)
{
  int t, r, f1, e, a, l, d;

  for(t = 0; t < MT; t++) {
    for(r = 0; r <=MR; r++) {
      for(f1 = 0; f1 <= MF; f1++) {
	for(e = 0; e <= ME; e++) {
	  for(a = 0; a <= MA; a++) {
	    for(l = 0; l < ML; l++) {
	      for(d = 0; d < MDEC; d++) {
		o_uvalue[t][r][f1][e][a][l][d] = 
		  uvalue[t][r][f1][e][a][l][d] = 
		  (1.0-P[0].simit) * o_uvalue[t][r][f1][e][a][l][d] +
		  P[0].simit * uvalue[t][r][f1][e][a][l][d];
		o_action[t][r][f1][e][a][l][d] = 
		  action[t][r][f1][e][a][l][d] = 
		  (1.0-P[0].simit) * o_action[t][r][f1][e][a][l][d] +
		  P[0].simit * action[t][r][f1][e][a][l][d];
	      }
		o_moult1[t][r][f1][e][a][l] = 
		  moult1[t][r][f1][e][a][l] = 
		  (1.0-P[0].simit) * o_moult1[t][r][f1][e][a][l] +
		  P[0].simit * moult1[t][r][f1][e][a][l];
	    }
	  }
	}
      }
    }
  }

}
/************************* copy_policy *************************/

