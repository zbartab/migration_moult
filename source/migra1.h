
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<fpu_control.h>
#include<time.h>

#define YES 1
#define NO 0



/*states*/
#define MR 12			/* max reserves */
#define MF 23			/* max feather */
#define MT 52			/* max week */
#define ME 2			/* max experience */
#define MA 8			/* max brood age */
#define MOL 13			/* moult length */
#define ML 2			/* max location */
#define MDEC 4			/* max decisions */

#define N_PATCH 5

#define TOPR ((double) MR)
#define TOPF ((double) MF-MOL)

#define UGRAIN 10
#define UTOL 1.0e-3

#define PI 3.1415926535898

#define MAX 900
#define SUM 901
#define ALL 902
#define BACKWARD 903
#define FORWARD 904
#define YOUNGS 905
#define NOYOUNGS 906
#define INIT 907
#define COPY 908
#define NORM 909
#define NONORM 910

/*actions*/
#define DELAY 1
#define START 2
#define KEEP 3
#define ABANDON 4
#define ABORT 5
#define FLYNORTH 6
#define FLYSOUTH 7
#define NOMOULT 8
#define MOULT1 9
#define BIRTOK 10

/* decision codes */
/*  if (a==0) */
#define D_DELAY 0
#define D_BIRTOK 1
#define D_FLYNORTH 2
#define D_FLYSOUTH 3

/*  if (a==1) */
#define D_DELAY 0
#define D_START 1
#define D_FLYNORTH 2
#define D_FLYSOUTH 3

/* if (a>1 && a<MA) */
#define D_ABORT 0
#define D_KEEP 1

/* if (a==MA) */
#define D_ABORT 0
#define D_ABANDON 1




typedef struct {	/***** feather stuff *****/
  int fixed;		/* abrasion fixed? 1: yes */
  double fora;		/* effect of foraging intensity on feather abrasion */
  double foraM;		/* effect of foraging intensity on moult speed */
  double flight;	/* effect of feather quaility on the flight ability */
  double worn;		/* the flight ability of birds with feathers in very poor condition */
  double NU;		/* stochasticity of moult length  */
} t_f;

typedef struct {	/***** migration stuff *****/
  double feather;	/* effect of migration on feather quality */
  double featherM;	/* effect of migration on moulting feather quality */
  double energy;	/* energetic cost of migration with top quality feathers */
  double preda;		/* predation risk during migration */
} t_m;

typedef struct {	/***** energetic stuff *****/
  double fora;		/* cost of foraging */
  double propF;		/* the proportion of flight used during foraging */
  double propNF;	/* 1-propF; just to save computation time */
  double def;		/* the cost of decreased flight eff. during foraging */
  double start;		/* cost of starting a brood */
  double keep;		/* amount of food needed by the brood */
  double birtok;	/* cost of territory occupation */
  double moult1;	/* cost of moulting feather 1 alone */
  double mass;		/* mass dependent cost */
  double baseM;		/* mean maintance cost */
  double baseS;		/* seasonality of maintance cost */
  double base[MT+1];	/* time dependent maintance cost */
} t_c;

typedef struct {	/***** predation stuff *****/
  double base;		/* mass independent predation risk */
  double mass;		/* mass dependent predation risk */
} t_p;

typedef struct {	/***** land stuff *****/
  double THETA;		/* parameter describing young foraging efficiency */
  double FOOD;		/* A_{food}: average amount of food over the year */
  double MAXEPS;	/* \epsilon: seasonality in food */
  double N0;		/* population size  */
  double DD;		/* measure of the strength of density dependence */
  double p;		/* probability parameter of the binom food distribution */
  double omega;		/* scaling factor for the food distribution */
  double pe;		/* probability of changing experience class */
  double n;		/* brood size */
  double n_birtok;	/* total number of territories */
  double K;		/* parameter of territory searching */
  double al;		/* parameter of territory searching */
} t_l;

typedef struct {
  t_f F;
  t_c C;
  t_p P;
  t_m M;
  t_l L;
  double tol;
  double alpha;
  double alpha2;
  double n_alpha;
  double K;
  double simit;
} type_par;


typedef struct {
  int r;
  int f1;
  int e;
  int a;
  int l;
  int alive;
  int pred;
  int birth;
} t_ind;
