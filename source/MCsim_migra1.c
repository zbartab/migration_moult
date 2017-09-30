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
#include"MCsim_migra1.h"

#define NBIRDS 2000
#define NYEARS 21
#define FNBIRD ((float) NBIRDS)

extern double ******states;
extern double uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
extern double moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];
extern double theta[ME+1][ML];
extern double FOOD[MT+1][ML];
extern double G[MT+1][ME+1][ML][N_PATCH];
extern double p_birtok[MT+1][ML];
extern double p_i[ML][N_PATCH];

extern type_par P[ML];

extern FILE *OutPut3;

double oFOOD[MT+1][ML];
type_par oP[ML];

long ID_base;

extern float randf(float);
extern int randi(int);
extern double newres(double *, int, int, int, int, int, double, int);
extern double newfq1(double, int, int);
extern void Migration(double *, double *, double *, int, int, int);
extern double P_birtok(double, double, int);
/*extern void Treatment(int, t_ind *);*/
extern void Food_Treatment(int, int);
extern double food_dens(double n, int l);

void ynormalise(double *, int );
void simulation(char *filename, int tt, t_manip *f_manip);
void MoultStr(int, char *);
void ActionStr(int migact, int act, char *S);
void WriteOut(int WRITE, int year, int t, char *S0, char *S1, char *S3,
	      t_bird *Bird, double u, int raf, int f1af, int laf, FILE *OutPut);
int smooth_r(double);
int smooth_f1(double);
int my_rdev(double *, int);
void free_bird(t_bird *bird);
t_bird *init_birds(int t);
t_bird *new_bird(t_bird *prev_bird);
t_bird *remove_bird(t_bird *bird);
void new_chick(t_bird *parent);

t_manip *ReadManip(char *filename);
void Manipulate(int year, int week, t_bird *f_bird, t_manip *f_manip);
void free_manip(t_manip *manip);
void Copy_Food(void);

int comp_week(int w1, int w2);
void set_state(t_bird *b, t_action ID, char a, int week, double u, 
	       int r_af, int f1_af);
void Record(int year, int t, char *S0, char *S1, char *S3,
	    t_bird *Bird, double u, int raf, int f1af);
void start_year(t_bird *f_bird);
void end_year(t_bird *f_bird, int year, char *filename);


/************************* start_year *************************/
void start_year(t_bird *f_bird)
{
  t_bird *bird;
  int l, a;

  for(bird = f_bird->next; bird != NULL; bird = bird->next){
    if(bird->e==ME) bird->record = YES;
    else bird->record = NO;
    for(l=0; l < ML; l++) {
      for(a=0; a < L_ACTIONS; a++) {
	bird->stat[l][a].count = 0;
	bird->stat[l][a].week = my_NA;
	bird->stat[l][a].r = my_NA;
	bird->stat[l][a].f1 = my_NA;
	bird->stat[l][a].r_af = my_NA;
	bird->stat[l][a].f1_af = my_NA;
	bird->stat[l][a].u = (double) my_NA;
	bird->stat[l][a].act = 'X';
      }
      bird->stat[l][DEATH].act = 'A';
    }
  }
}
/************************* start_year *************************/

/************************* end_year *************************/
void end_year(t_bird *f_bird, int year, char *filename)
{
 t_bird *bird;
 int l, a, i, h;
 double hd;
 char filen[256];
 FILE *f_apd;
 char nev[30];

  sprintf(filen,"%s-%d.apd",filename,year);
  if ((f_apd = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }

  for(a=0; a < L_ACTIONS; a++) {
    for(l=0; l<ML; l++) {
      switch (a) {
      case DEATH:
	if(l==0) fprintf(f_apd,"d.act d.loc d.n d.w d.res d.fea d.u ");
	break;
      case WINTER:
	if(l==0) fprintf(f_apd,"w.loc w.n w.w w.res w.fea w.u ");
	break;
      case E_WINTER:
	if(l==0) fprintf(f_apd,"l.loc l.n l.w l.res l.fea l.u ");
	break;
      case SUMMER:
	if(l==0) fprintf(f_apd,"u.loc u.n u.w u.res u.fea u.u ");
	break;
      case TERR:
	fprintf(f_apd,"n.te%d w.te%d r.te%d f.te%d u.te%d ",
		l,l,l,l,l);
	break;
      case BREED:
	fprintf(f_apd,"n.br%d w.br%d r.br%d f.br%d u.br%d ",
		l,l,l,l,l);
	break;
      case SU_BR:
	fprintf(f_apd,"n.su%d w.su%d r.su%d f.su%d u.su%d ",
		l,l,l,l,l);
	break;
      case MOULT:
	fprintf(f_apd,"n.mo%d w.mo%d r.mo%d f.mo%d u.mo%d ",
		l,l,l,l,l);
	break;
      case E_MOULT:
	fprintf(f_apd,"n.en%d w.en%d r.en%d f.en%d u.en%d ",
		l,l,l,l,l);
	break;
      case NORTH:
	fprintf(f_apd,"n.no%d w.no%d r.no%d f.no%d u.no%d ar.no%d af.no%d ",
		l,l,l,l,l,l,l);
	break;
      case SOUTH:
	fprintf(f_apd,"n.so%d w.so%d r.so%d f.so%d u.so%d ar.so%d af.so%d ",
		l,l,l,l,l,l,l);
	break;
      }
    }
  }
  fprintf(f_apd,"\n");

 for(bird = f_bird->next; bird != NULL; bird = bird->next){
    if(bird->record==YES){
      for(a = 0; a < L_ACTIONS; a++) {
	if(a==DEATH) {
	  i = -1;
	  for(l = 0; l < ML; l++) {
	    if(bird->stat[l][a].act != 'A') i = l;
	  }
	  l = (i == -1 ? 0 : i);
	  fprintf(f_apd,"%c ", bird->stat[l][a].act);
	  fprintf(f_apd,"%d ", l);
	  if((h=bird->stat[l][a].count)!=my_NA) fprintf(f_apd,"%d ",h);
	  else fprintf(f_apd,"NA ");
	  if((h=bird->stat[l][a].week)!=my_NA) fprintf(f_apd,"%d ",h);
	  else fprintf(f_apd,"NA ");
	  if((h=bird->stat[l][a].r)!=my_NA) fprintf(f_apd,"%d ",h);
	  else fprintf(f_apd,"NA ");
	  if((h=bird->stat[l][a].f1)!=my_NA) fprintf(f_apd,"%d ",h);
	  else fprintf(f_apd,"NA ");
	  if((hd=bird->stat[l][a].u)!=my_NA) fprintf(f_apd,"%lf ",hd);
	  else fprintf(f_apd,"NA ");
	  //	  fprintf(f_apd,"%d %d %d %d %lf ", bird->stat[l][a].count,
	  //  bird->stat[l][a].week, bird->stat[l][a].r,
	  //  bird->stat[l][a].f1, bird->stat[l][a].u);
	} else {
	  i = -1;
	  for(l = 0; l < ML; l++) {
	    if (a==WINTER || a==E_WINTER || a==SUMMER) {
	      if(bird->stat[l][a].count > 0) {
		fprintf(f_apd,"%d ", l);
		if((h=bird->stat[l][a].count)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		if((h=bird->stat[l][a].week)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		if((h=bird->stat[l][a].r)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		if((h=bird->stat[l][a].f1)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		if((hd=bird->stat[l][a].u)!=my_NA) fprintf(f_apd,"%lf ",hd);
		else fprintf(f_apd,"NA ");
		//fprintf(f_apd,"%d %d %d %d %lf ", bird->stat[l][a].count,
		//bird->stat[l][a].week, bird->stat[l][a].r,
		//bird->stat[l][a].f1, bird->stat[l][a].u);
		i = 1;
	      } else if(l == ML-1 && i == -1) {
		fprintf(f_apd,"NA NA NA NA NA NA ");
		//fprintf(f_apd,"%d ", my_NA);
		//fprintf(f_apd,"%d %d %d %d %lf ", bird->stat[0][a].count,
		//bird->stat[0][a].week, bird->stat[0][a].r,
		//bird->stat[0][a].f1, bird->stat[0][a].u);
	      }
	    } else {
	      if((h=bird->stat[l][a].count)!=my_NA) fprintf(f_apd,"%d ",h);
	      else fprintf(f_apd,"NA ");
	      if((h=bird->stat[l][a].week)!=my_NA) fprintf(f_apd,"%d ",h);
	      else fprintf(f_apd,"NA ");
	      if((h=bird->stat[l][a].r)!=my_NA) fprintf(f_apd,"%d ",h);
	      else fprintf(f_apd,"NA ");
	      if((h=bird->stat[l][a].f1)!=my_NA) fprintf(f_apd,"%d ",h);
	      else fprintf(f_apd,"NA ");
	      if((hd=bird->stat[l][a].u)!=my_NA) fprintf(f_apd,"%lf ",hd);
	      else fprintf(f_apd,"NA ");
	      //fprintf(f_apd,"%d %d %d %d %lf ", bird->stat[l][a].count,
	      //      bird->stat[l][a].week, bird->stat[l][a].r,
	      //      bird->stat[l][a].f1, bird->stat[l][a].u);
	      if(a==NORTH || a==SOUTH) {
		if((h=bird->stat[l][a].r_af)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		if((h=bird->stat[l][a].f1_af)!=my_NA) fprintf(f_apd,"%d ",h);
		else fprintf(f_apd,"NA ");
		//fprintf(f_apd,"%d %d ", bird->stat[l][a].r_af,
		//	bird->stat[l][a].f1_af);
	      }
	    }
	  }
	}
      }
      fprintf(f_apd,"\n");
    }
   if(bird->remove == YES) bird = remove_bird(bird);
 }
 fclose(f_apd); 
}
/************************* end_year *************************/

/************************* simulation *************************/
void simulation(char *filename, int tt, t_manip *f_manip)
{
  t_bird *bird;
  char filen[256], S1[30], S2[30], S3[30];
  int t, to, r, f1, e, a, act, migact, ido, l;
  double lambda, u, pmort, e_exp;
  int year, i;
  FILE *OutPut1;
  double rf, ff1, p0, p1, keep;
  int alive;
  int WRITE;
  double n_fora[ML], g, N_birtok[ML], N_bird[ML];
  t_bird *First;

  ID_base = randi(10000000);
  First = init_birds(tt);

  for(t = 0; t <= MT; t++)
    for(l = 0; l < ML; l++)
      oFOOD[t][l] = FOOD[t][l];
  for(l = 0; l < ML; l++) oP[l] = P[l];

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".dat");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  /*  fprintf(OutPut1, 
      "year\tweek\tbird\tevent\tmoult1\taction\tres\tfea1\texper\tbrood\tloc\tuval\n");*/

  for(year = 0; year < NYEARS; year++){
    alive = 0;
    WRITE = NO;//year == 20;// || year == 21 || year == 19;
    start_year(First);
    for(to = 0; to < MT; to++){
      t = (to+tt) % MT;
      ido = year*52 + t;
      /*      if (tt != -1) Food_Treatment(ido, t);*/
      Copy_Food();
      for(l = 0; l < ML; l++) {
	P[l] = oP[l];
	N_birtok[l] = N_bird[l] = n_fora[l] = 0.0;
      }
      Manipulate(year,t,First,f_manip);
      for(bird = First->next; bird != NULL; bird = bird->next){
	if(bird->remove==YES) continue;
	migact = 0;
	keep = 0.0;
	pmort = 0.0;
	if (bird->r <= 0.0) {
	  if (bird->mig_died==NO) {
	    strcpy(S2,"S");
	    WriteOut(WRITE, year,t,S2,S1,S3,bird,u,-9,-9,-9,OutPut1);
	  }
	  if(bird->record == YES) {
	    bird->remove = YES;
	    continue;
	  } else {
	    bird = remove_bird(bird);
	  }
	} else if (bird->pred) {
	  if (bird->mig_died==NO) {
	    strcpy(S2,"P");
	    WriteOut(WRITE,year,t,S2,S1,S3,bird,u,-9,-9,-9,OutPut1);
	  }
	  if(bird->record == YES) {
	    bird->remove = YES;
	    continue;
	  } else {
	    bird = remove_bird(bird);
	  }
	} else {
	  r = bird->r;
	  f1 = bird->f1;
	  e = bird->e;
	  a = bird->a;
	  l = bird->l;
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
	  if (randf(1.0) < pmort) bird->pred = 1;
	  if(r > 0 && bird->pred != 1) {
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
	    WriteOut(WRITE,year,t,"A",S1,S3,bird,u,r,f1,l,OutPut1);
	    n_fora[l] += u * theta[e][l];
	    if(a>0) N_birtok[l] += 1.0;
	    else if(act==BIRTOK) N_bird[l] += 1.0;
	  } else {
	    bird->mig_died = YES;
	    ActionStr(migact,DELAY,S3);
	    if (r <= 0.0) {
	      strcpy(S2,"SM");
	    } else if (bird->pred) {
	      strcpy(S2,"PM");
	    } else {
	      strcpy(S2,"X");
	    }
	    WriteOut(WRITE, year,t,S2,S1,S3,bird,u,r,f1,l,OutPut1);
	  }
	  bird->r = r;
	  bird->f1 = f1;
	  bird->e = e;
	  bird->a = a;
	  bird->l = l;
	  bird->act = act;
	  bird->u = u;
	}
      }
      for(l = 0; l < ML; l++) {
	p_birtok[t][l] = P_birtok(N_bird[l],N_birtok[l],l);
	for(e = 0; e <= ME; e++) {
  	  g = (FOOD[t][l]*food_dens(n_fora[l],l)) * theta[e][l];  
	  for(i = 0; i < N_PATCH; i++) {
	    G[t][e][l][i] =  TOPR * 
	      (g - P[l].L.p*(N_PATCH-1)*P[l].L.omega + i*P[l].L.omega);
	  }
	}
      }
      for(bird = First->next; bird != NULL; bird = bird->next){
	if(bird->remove==YES) continue;
	pmort = 0.0;
	if (t==0) alive++;
	if (bird->r > 0.0 && bird->pred != 1) {
	  r = bird->r;
	  f1 = bird->f1;
	  e = bird->e;
	  a = bird->a;
	  l = bird->l;
	  act = bird->act;
	  u = bird->u;
	  if(act == KEEP || act == ABANDON) keep = P[l].C.keep;
	  else keep = 0.0;
	  e_exp = newres(&pmort,t,r,f1,e,l,u,act);
	  ff1 = newfq1(e_exp,f1,l);
	  rf = ((double) r) + u*G[t][e][l][my_rdev(p_i[l],N_PATCH)] - 
	    e_exp - keep; 
	  bird->r = smooth_r(rf);
	  bird->f1 = smooth_f1(ff1);
	  if (randf(1.0) < pmort) bird->pred = 1;
	  if (randf(1.0)<P[l].L.pe) bird->e = (e+1 < ME ? e+1 : ME);
	  bird->a = a;
	  if (act == START) bird->a = 2;
	  else if (act == KEEP)  bird->a = (a+1 > MA ? MA : a+1);
	  else if (act == BIRTOK && randf(1.0) < p_birtok[t][l]) 
	    bird->a = 1;
	  else if (act == ABORT) bird->a = 1;
	  else if (act == ABANDON) {
	    bird->a = 1;
	    for(i = 0; i < P[l].L.n; i++) {
	      new_chick(bird);
	    }
	  }
	  bird->l = l;
	}
      }
    }
    printf("%d, %d\n", year, alive);
    fprintf(OutPut3,"%d, %d\n", year, alive);
    end_year(First,year,filename);
  }
  fclose(OutPut1);
}
/************************* simulation *************************/

/************************* new_chick *************************/
void new_chick(t_bird *parent)
{
  static long id = 0;
  t_bird *chick;

  if (id==0) id = ID_base + NBIRDS + 1;

  chick = new_bird(parent);
  chick->ID = id++;
  chick->alive = 1;
  chick->pred = 0;
  chick->mig_died = NO;
  chick->r = 0.5*TOPR;
  chick->f1 = MF;
  chick->e = 0;
  chick->a = 0;
  chick->act = DELAY;
  chick->u = 0.0;
  chick->l = parent->l;
}
/************************* new_chick *************************/

/************************* ynormalise *************************/
t_bird *init_birds(int t)
{
  int r, f1, e, a, l;
  double lambda, tmp, L1, L2;
  int i;
  t_bird *bird, *next_bird, *f_bird;

  f_bird = bird = new_bird(NULL);
  lambda = 0.0;
  for(r = 1; r <= MR; r++)
  for(f1 = 0; f1 <= MF; f1++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++)
  for(l = 0; l < ML; l++){
    lambda += states[t][r][f1][e][a][l];
  }
  L1 = L2 = 0.0;
  i = 0;
  if (lambda > 0.0){
    for(r = 1; r <= MR; r++)
    for(f1 = 0; f1 <= MF; f1++)
    for(e = 0; e <= ME; e++)
    for(a = 0; a <= MA; a++)
    for(l = 0; l < ML; l++){
      L2 += (states[t][r][f1][e][a][l] /= lambda);
      while(((((float) i)/FNBIRD) < L2) && (i<NBIRDS)){
	/*	next_bird = new_bird(bird);
		bird->next = next_bird;*/
	bird = new_bird(bird);
	/*	bird = next_bird;*/
	bird->ID = ID_base + i;
	bird->alive = 1;
	bird->pred = 0;
	bird->mig_died = NO;
	bird->r = r;
	bird->f1 = f1;
	bird->e = e;
	bird->a = a;
	bird->l = l;
	i++;
      }
    }
  }
  else {
    fprintf(OutPut3,"ERROR: Something is wrong with lambda in `init_bird'\n");
    exit(1);
  }
  return f_bird;
}
/************************* end of ynormalise ******************/

/************************* new_bird *************************/
t_bird *new_bird(t_bird *prev_bird)
{
  t_bird *new_b, *next_bird;
  int l, a;
  
  new_b = (t_bird *) malloc(sizeof(t_bird));
  if(new_b == NULL) {
    fprintf(OutPut3,"ERROR: not enough memory in `new_bird'\n");
    exit(1);
  }
  next_bird = NULL;
  if(prev_bird != NULL) {
    next_bird = prev_bird->next;
    prev_bird->next = new_b;
  }
  new_b->prev = prev_bird;
  if(next_bird != NULL) {
    next_bird->prev = new_b;
  }
  new_b->next = next_bird;
  new_b->remove = NO;
  new_b->record = NO;
  for(l=0; l < ML; l++) {
    for(a=0; a < L_ACTIONS; a++) {
      new_b->stat[l][a].count = 0;
      new_b->stat[l][a].week = my_NA;
      new_b->stat[l][a].r = my_NA;
      new_b->stat[l][a].f1 = my_NA;
      new_b->stat[l][a].r_af = my_NA;
      new_b->stat[l][a].f1_af = my_NA;
      new_b->stat[l][a].u = (double) my_NA;
      new_b->stat[l][a].act = 'X';
    }
  }
  return new_b;
}
/************************* new_bird *************************/

/************************* remove_bird *************************/
t_bird *remove_bird(t_bird *bird)
{
  t_bird *prev, *next;

  if(bird==NULL) return;
  prev = bird->prev;
  next = bird->next;
  
  if(next != NULL) {
    next->prev = prev;
    bird->next = NULL;
  }
  prev->next = next;
  free_bird(bird);
  return prev;
}
/************************* remove_bird *************************/

/************************* free_bird *************************/
void free_bird(t_bird *bird)
{
  if(bird == NULL) return;
  if(bird->next != NULL) free_bird(bird->next);
  free(bird);
  bird = NULL;
  return;
}
/************************* free_bird *************************/

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
void WriteOut(int WRITE, int year, int t, char *S0, char *S1, char *S3,
	      t_bird *Bird, double u, int raf, int f1af, int laf, FILE *OutPut)
     /* this function writes out the bird's state and its decision at
	the begining of week t; .af variables are after migration variables*/
{
  if( (t==0) || (t==26) || (t==51) || 
      !(strcmp(S0,"A")==0 && strcmp(S1,"N")==0 && strcmp(S3,"D")==0) ) {
    if (WRITE) 
      fprintf(OutPut,
	      "%d\t%d\t%X\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\n",
	      year, t, Bird->ID, S0, S1, S3, Bird->r, Bird->f1, Bird->e,Bird->a,
	      Bird->l,u,raf,f1af,laf);
    Record(year, t, S0, S1, S3, Bird, u, raf, f1af);
  }
  
}
/************************* WriteOut *************************/

/************************* comp_week *************************/
int comp_week(int w1, int w2)
{
  if (abs(w1-w2) < 26) 
    return w1 < w2;
  else
    return w1 > w2;
}
/************************* comp_week *************************/


/************************* set_state *************************/
void set_state(t_bird *b, t_action ID, char a, int week, double u, 
	       int r_af, int f1_af)
{
  int l;

  l = b->l;
  if(b->stat[l][ID].count==0) {
    b->stat[l][ID].count = 1;
    b->stat[l][ID].act = a;
    b->stat[l][ID].week = week;
    b->stat[l][ID].r = b->r;
    b->stat[l][ID].f1 = b->f1;
    b->stat[l][ID].r_af = r_af;
    b->stat[l][ID].f1_af = f1_af;
    b->stat[l][ID].u = u;
  } else if (comp_week(week,b->stat[l][ID].week)) {
    b->stat[l][ID].count += 1;
    b->stat[l][ID].act = a;
    b->stat[l][ID].week = week;
    b->stat[l][ID].r = b->r;
    b->stat[l][ID].f1 = b->f1;
    b->stat[l][ID].r_af = r_af;
    b->stat[l][ID].f1_af = f1_af;
    b->stat[l][ID].u = u;
  } else b->stat[l][ID].count += 1;
}
/************************* set_state *************************/


/************************* Record *************************/
void Record(int year, int t, char *S0, char *S1, char *S3,
	      t_bird *Bird, double u, int raf, int f1af)
{
  if(t==0) set_state(Bird, WINTER, 'W', t, u, my_NA, my_NA);
  if(t==26) set_state(Bird, SUMMER, 'U', t, u, my_NA, my_NA);
  if(strcmp(S0,"P")==0) {
    set_state(Bird, DEATH, 'P', t, u, my_NA, my_NA);
  } else if(strcmp(S0,"S")==0) {
    set_state(Bird, DEATH, 'F', t, u, my_NA, my_NA);
  } else if(strcmp(S0,"PM")==0) {
    set_state(Bird, DEATH, 'R', t, u, my_NA, my_NA);
  } else if(strcmp(S0,"SM")==0) {
    set_state(Bird, DEATH, 'H', t, u, my_NA, my_NA);
  } else {
    if(strcmp(S1,"S")==0) {
      set_state(Bird, MOULT, 'M', t, u, my_NA, my_NA);
    } else if(strcmp(S1,"E")==0) {
      set_state(Bird, E_MOULT, 'e', t, u, my_NA, my_NA);
    }
    if(strstr(S3,"No")==S3) {
      set_state(Bird,NORTH,'N',t, u,raf,f1af);
      S3 += 2;
    }
    if(strstr(S3,"So")==S3) {
      set_state(Bird,SOUTH,'S',t, u,raf,f1af);
      S3 += 2;
    }
    Bird->r = raf;
    Bird->f1 = f1af;
    if(strcmp(S3,"S")==0) 
      set_state(Bird,BREED,'B',t, u,my_NA,my_NA);
    else if(strcmp(S3,"A")==0) 
      set_state(Bird,SU_BR,'a',t, u,my_NA,my_NA);
    else if(strcmp(S3,"B")==0) 
      set_state(Bird,TERR,'T',t, u,my_NA,my_NA);
  }
  if(t==51) set_state(Bird, E_WINTER, 'L', t, u, my_NA, my_NA);
}
/************************* Record *************************/

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

/************************* new_manip *************************/
t_manip *new_manip(t_manip *prev_manip)
{
  t_manip *new_b, *next_manip;
  
  new_b = (t_manip *) malloc(sizeof(t_manip));
  if(new_b == NULL) {
    fprintf(OutPut3,"ERROR: not enough memory in `new_manip'\n");
    exit(1);
  }
  next_manip = NULL;
  if(prev_manip != NULL) {
    next_manip = prev_manip->next;
    prev_manip->next = new_b;
  }
  new_b->prev = prev_manip;
  if(next_manip != NULL) {
    next_manip->prev = new_b;
  }
  new_b->next = next_manip;
  return new_b;
}
/************************* new_manip *************************/

/************************* remove_manip *************************/
t_manip *remove_manip(t_manip *manip)
{
  t_manip *prev, *next;

  if(manip==NULL) return;
  prev = manip->prev;
  next = manip->next;
  
  if(next != NULL) {
    next->prev = prev;
    manip->next = NULL;
  }
  prev->next = next;
  free_manip(manip);
  return prev;
}
/************************* remove_manip *************************/

/************************* free_manip *************************/
void free_manip(t_manip *manip)
{
  if(manip == NULL) return;
  if(manip->next != NULL) free_manip(manip->next);
  free(manip);
  manip = NULL;
  return;
}
/************************* free_manip *************************/

/************************* ReadManip *************************/
t_manip *ReadManip(char *filename)
{
  FILE *InPut;
  char filen[30];
  char *ext = ".man";
  int errors = 0;
  char junk[50];
  int m_type, n_manip;
  int m_year, m_week, m_loc;
  double m_vol;
  t_manip *manip, *f_manip;

  f_manip = manip = new_manip(NULL);
  strcpy(filen,""); strcat(filen,filename); strcat(filen,ext);
  if ((InPut = fopen(filen,"rt")) == NULL) {
    fprintf(OutPut3,"WARNING: Can't open the manip file: %s\n", filen);
    /*    exit (1);*/
    return f_manip;
  }
  errors = 0;
  n_manip = 0;
  while(fscanf(InPut,"%s%d%d%d%lf",
	       junk,&m_year,&m_week,&m_loc,&m_vol) != EOF) {
    n_manip++;
    if(strcmp(junk,"preda")==0) m_type = M_PREDA;
    else if(strcmp(junk,"food")==0) m_type = M_FOOD;
    else if(strcmp(junk,"reserve")==0) m_type = M_RESERVE;
    else if(strcmp(junk,"feather")==0) m_type = M_FEATHER;
    else if(strcmp(junk,"remove")==0) m_type = M_BROOD_REMOVE;
    else if(strcmp(junk,"elong")==0) m_type = M_BROOD_ELONG;
    else if(strcmp(junk,"terr")==0) m_type = M_TERR;
    else errors++;
    if(m_year < 0) errors++;
    if(m_week < 0 || m_week > 51) errors++;
    if(m_loc <0 || m_loc >= ML) errors++;
    if(errors > 0) {
      fprintf(OutPut3,"WARNING: Manipulation %d is invalid\n", n_manip);
    } else {
      manip = new_manip(manip);
      manip->type = m_type;
      manip->year = m_year;
      manip->week = m_week;
      manip->loc = m_loc;
      manip->vol = m_vol;
    }
  }
  fclose(InPut);
  return f_manip;
}
/************************* ReadManip *************************/

/************************* Copy_Food *************************/
void Copy_Food(void)
{
  int t, l;
  
  for(t = 0; t <= MT; t++) 
    for(l = 0; l < ML; l++) {
      FOOD[t][l] = oFOOD[t][l];
    }
}
/************************* Copy_Food *************************/

/************************* Manipulate *************************/
void Manipulate(int year, int week, t_bird *f_bird, t_manip *f_manip)
{
  t_manip *manip;
  t_bird *bird;

  double f;
  
  for(manip = f_manip->next; manip != NULL; manip = manip->next) {
    if(manip->year == year && manip->week == week) {
      switch(manip->type) {
      case M_PREDA:
	P[manip->loc].P.mass *= manip->vol;
	if(P[manip->loc].P.mass < 0.0) P[manip->loc].P.mass = 0.0;
	break;
      case M_FOOD:
	FOOD[week][manip->loc] *= manip->vol;
	if(FOOD[week][manip->loc] < 0.0) FOOD[week][manip->loc] = 0.0;
	break;
      case M_RESERVE:
	for(bird = f_bird->next; bird != NULL; bird = bird->next){
	  if (bird->l == manip->loc) {
	    bird->r = (int) (((double)bird->r) + (manip->vol)*TOPR);
	  }
	}
	break;
      case M_FEATHER:
	for(bird = f_bird->next; bird != NULL; bird = bird->next){
	  if (bird->l == manip->loc) {
	    bird->f1 = (int) (((double)bird->f1) + (manip->vol)*TOPF);
	    if(bird->f1 < 0) bird->f1 = 0;
	  }
	}
	break;
      case M_BROOD_ELONG:
      case M_BROOD_REMOVE:
      case M_TERR:
	break;
      }
    } else if(manip->year < year) manip = remove_manip(manip);
  }
}
/************************* Manipulate *************************/



