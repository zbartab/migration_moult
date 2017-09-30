

#include"migra1.h"


double repval[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double youngs[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double uvalue[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    action[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult1[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult2[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

type_par P;


double pa = 1.0; /* factor affecting the length of a breeding attempt */

long idum=(-13);

FILE *OutPut3;

typedef struct {
  int tt;
  double rr;
  double ff;
} t_tr;

t_tr Treats;

extern void Ind(char *, int);
extern void iniFP(void);
extern void inigain(double);

void ReadIn(char *);
void BinRead(char *);
void hiba(char *);
void Treatment(double *, double *);

/************************* main *************************/
int main(int argc, char *argv[])
{
  char filename[30];
  char filen[30];
  int i;
  float differ, pdiff;
  float In, D;

  Treats.tt = -1;
  Treats.rr = 0.0;
  Treats.ff = 0.0;
  idum = (long) time(NULL);
  if (idum > 0) idum = -idum;

  if ( argc > 4 )
    {
      strcpy(filename, argv[1]);
    }
  else strcpy(filename, "moult");
  if(argc > 3) {
    Treats.tt = atoi(argv[argc-3]);
    Treats.rr = atof(argv[argc-2]);
    Treats.ff = atof(argv[argc-1]);
  }
  else {
    hiba("Usage: moult_treat [database] tt rr ff\nwhere\ttt: the time of treatment\n\trr: increase of reserves\n\tff: feather manipulation\n");  
  }
  ReadIn(filename);
  inigain(1.0);
  iniFP();
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".sim");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    printf("Can't open the output file\n");
    exit (1);
  }
  fprintf(OutPut3,"Seed: %ld\n", idum);
  BinRead(filename);
  Ind(filename, Treats.tt);
  fclose(OutPut3);
}
/************************* end of main ******************/


/************************* BinRead *************************/
void BinRead(char *filename)
{
  int t, r, f1, e, a, l;
  char filen[30];
  FILE *BinOut;
  int n_y, n_a;
  int wrote;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".bin");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("Can't open the input file\n");
    exit (1);
  }
  if (fread(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  if (n_y != (MT+1)*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1)*LOCS) 
    hiba("Size of youngs is not correct");
  if (fread(repval, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2a");
  if (fread(youngs, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2b");
  if (fread(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinWrite 3");
  if (n_a != MT*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1)*LOCS)
    hiba("Size of action is not correct");
  if (fread(uvalue, sizeof(double),n_a,BinOut) != n_a) hiba("BinWrite 4");
  if (fread(action, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 5");
  if (fread(moult1, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 6");
  if (fread(moult2, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 7");
  fclose(BinOut);
}
/************************* BinRead *************************/

/************************* hiba *************************/
void hiba(char *leiras)
{
  fprintf(stderr, "%s\n", leiras);
  exit(1);
}
/************************* hiba *************************/


/************************* Treatment *************************/
void Treatment(double *rf, double *ff1)
{
  double s;

  *rf += Treats.rr;
  if (*ff1>=0.0) *ff1 = ((s=*ff1+Treats.ff)<0 ? 0.0 : s);
}
/************************* Treatment *************************/
