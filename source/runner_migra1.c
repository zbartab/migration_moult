
/*
  runner_migra1: individual based monte carlo simulation of migratory
  birds. It is based on the optimal annual routine model for a migratory
  bird (migra1)

  developers: Z. Barta       University of Bristol
              R.K. Welham    University of Bristol
	      J.M. McNamara  University of Bristol
	      A.I. Houston   University of Bristol

  maintainer: Z. Barta <zbarta@delfin.unideb.hu>

*/

#include"migra1.h"
#include"MCsim_migra1.h"

#define MAXN 100
#define MAX_BREED 10

double ucritcl[MT][ME+1][ML];
double unocare[MR+1][MF+1][ME+1][ML];
double vnocare[MR+1][MF+1][ME+1][ML];

double ******states;

double *****h_rv;		/* temporary storage because of migration */


double uvalue[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double action[MT][MR+1][MF+1][ME+1][MA+1][ML][MDEC];
double moult1[MT][MR+1][MF+1][ME+1][MA+1][ML];

extern double G[MT+1][ME+1][ML][N_PATCH];
extern double oG[MT+1][ME+1][ML][N_PATCH];
extern double theta[ME+1][ML];
extern double FOOD[MT+1][ML];

type_par P[ML];

long idum=(-13);

int Breed[MAX_BREED][2][ML];
int n_t = 0;

FILE *OutPut3;

extern void iniFP(void);
extern void inigain(double);
extern void simulation(char *, int, t_manip *f_manip);
extern void ReadIn(char *);
extern void iniBinom(void);
extern void free_arrays(double ******rv);
extern double ******init_arrays(void);
extern void free_h_rv(double *****rv);
extern double *****init_h_rv(void);
extern void hiba(char *);
extern t_manip *ReadManip(char *filename);
extern void free_manip(t_manip *manip);

void BinRead(char *filename);
double get_treat(int argc, char *argv[], char *filename, char *ext, int *n);
void G_out(char *filename);

/************************* main *************************/
int main(int argc, char *argv[])
{
  char filename[128];
  char filen[128];
  char ext[128];
  double lambda_back, diff_back;
  double lambda_for = 1.0, diff_for;
  double prev3lam;
  int I, l, i;
  int Kapcsolo = 0;
  double simit;
  double n;
  t_manip *f_manip;
  double P_size;

  idum = (long) time(NULL);
  if (idum > 0) idum = -idum;

  states = init_arrays();
  h_rv = init_h_rv();

  P_size = get_treat(argc, argv, filename, ext, &n_t);
  strcpy(filen,""); strcat(filen,filename); strcat(filen,"-"); 
  strcat(filen,ext);
  strcat(filen,".sum");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    fprintf(OutPut3,"ERROR: Can't open the output file\n");
    exit (1);
  }
  fprintf(OutPut3,"filename: %s\n", filename);
  fprintf(OutPut3,"Seed: %ld\n", idum);
  ReadIn(filename);
  for(l=0; l < ML; l++) P[l].L.N0 *= P_size;
  iniFP();
  iniBinom();
  inigain(1.0);
  BinRead(filename);
  strcat(filename,"-");strcat(filename,ext);
  f_manip = ReadManip(filename);
  simulation(filename, 0, f_manip);
  G_out(filename);
  free_manip(f_manip);
  free_h_rv(h_rv);
  free_arrays(states);
  fclose(OutPut3);
  return 0;
}
/************************* end of main ******************/


/************************* BinRead *************************/
void BinRead(char *filename)
{
  char filen[30];
  FILE *BinOut;
  int n_y, n_a, n_m;
  int t, r, f1, e, a, l;

  n_y = (MT+1)*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML;
  n_m = MT*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML;
  n_a = MT*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML*MDEC;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".food");
  if ((BinOut = fopen(filen,"r")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  for(t = 0; t <= MT; t++) {
    for(e = 0; e <= ME; e++) {
      for(l = 0; l < ML; l++) {
	fscanf(BinOut,"%f",&(G[t][e][l][(N_PATCH-1)/2]));
      }
    }
  }
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".stb");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fread(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinRead 1");
  if(n_y != (MT+1)*(MR+1)*(MF+1)*(ME+1)*(MA+1)*ML) {
    hiba("ERROR: size of states does not match");
  }
  for(t = 0; t <= MT; t++) {
  /*  t = 0;*/
  for(r = 0; r <=MR; r++) {
    for(f1 = 0; f1 <= MF; f1++) {
      for(e = 0; e <= ME; e++) {
	for(a = 0; a <= MA; a++) {
	  for(l = 0; l < ML; l++) {
	    if (fread(&(states[t][r][f1][e][a][l]), 
		      sizeof(double),1,BinOut) 
		!= 1) 
	      hiba("BinRead 2b");
	  }
	}
      }
    }
  }
  }
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".uvb");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fread(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinRead 3");
  if (fread(uvalue, sizeof(double),n_a,BinOut) != n_a) hiba("BinRead 4");
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".acb");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fread(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinRead 1");
  if (fread(action, sizeof(double),n_a,BinOut) != n_a) hiba("BinRead 5");
  fclose(BinOut);

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".mob");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  if (fread(&n_m,sizeof(int),1,BinOut) != 1) hiba("BinRead 1");
  if (fread(moult1, sizeof(double),n_m,BinOut) != n_m) hiba("BinRead 6");
  fclose(BinOut);
}
/************************* BinRead *************************/

/************************* get_treat *************************/
double get_treat(int argc, char *argv[], char *filename, char *ext, int *n) 
     /* args format: -b location start end  filename extension 
                     -N scalefactor_for_population_size */
{
  int l, i = 1, p = 0;
  int name_read = NO;
  double P;

  strcpy(filename,"baseline");
  strcpy(ext,"sim");
  P = 10.0;
  while (i < argc) {
    if (strcmp(argv[i],"-b") == 0) {
      l = atoi(argv[++i]); /* breeding location */
      Breed[p][0][l] = atoi(argv[++i]); /* breeding start */
      Breed[p][1][l] = atoi(argv[++i]); /* breeding finish */
      p++;
    } else if (strcmp(argv[i],"-N") == 0) {
      P = atoi(argv[++i]);
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
  return P;
}
/************************* end of get_treat ******************/


/************************* G_out *************************/
void G_out(char *filename) 
{
  FILE *OutPut;
  int t, e, l, n_p;
  char filen[500];

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".food");
  if ((OutPut = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  fprintf(OutPut,"week exper loc n.patch G\n");
  for(t = 0; t < MT; t++)
    for(e = 0; e <= ME; e++)
      for(l = 0; l < ML; l++)
	for(n_p = 0; n_p < N_PATCH; n_p++)
	  fprintf(OutPut, "%d %d %d %d %f\n", t, e, l, n_p, G[t][e][l][n_p]);
  fclose(OutPut);
}
/************************* G_out *************************/

