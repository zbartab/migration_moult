
#define my_NA -9


typedef struct BIRD t_bird;
typedef struct MANIP t_manip;
typedef struct STATES t_states;

typedef enum ACT t_action;

enum ACT {WINTER, SUMMER, E_WINTER, TERR, BREED, SU_BR, MOULT, 
	  E_MOULT, NORTH, SOUTH, DEATH};

#define L_ACTIONS 11		/* !!! don't forget to update if ACTIONS changed !!! */

struct STATES {
  char act;			/* the performed action */
  int count;			/* the count of action */
  int week;			/* the week of the first occurance */
  int r;			/* reserves when performed */
  int f1;			/* feather when performed */
  int r_af;			/* reserves after migration */
  int f1_af;			/* feather after migration */
  double u;			/* foraging intensity */
};

struct BIRD {		/* the bird data */
  t_bird *prev;
  t_bird *next;
  long ID;
  int alive;
  int r;
  int f1;
  int e;
  int a;
  int l;
  int pred;
  int act;
  double u;
  int mig_died;
  int remove;
  int record;
  t_states stat[ML][L_ACTIONS];
};

struct MANIP {
  t_manip *prev;
  t_manip *next;
  int type;
  int year;
  int week;
  int loc;
  double vol;
};


/* manipulation types */
#define M_INVALID 0		/* invalid type */
#define M_PREDA 1		/* predation risk */
#define M_FOOD 2		/* food supply */
#define M_RESERVE 3		/* reserves */
#define M_BROOD_REMOVE 4	/* brood removal a' = 0 */
#define M_BROOD_ELONG 5         /* brood elongation a' = a-1 */
#define M_TERR 6		/* territory */
#define M_FEATHER 7		/* feather quality */

#define M_YEAR_ZERO 0		/* manipulation only in zero year */
#define M_YEAR_REP 1		/* manipulation in every year after zero */
