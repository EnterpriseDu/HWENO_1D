/* 
 * This is a implementation of GRP scheme for 1-D
 * scalar conservation law:
 *                    u_t + f(u)_x = 0.
 *
 * 
 * The protential of reforming this implementation to solve
 * systems of equations of consevation laws is quite limited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "file_io.h"
#include "solver.h"
#include "reconstruction.h"
#include "Riemann_solver.h"




#ifndef N_CONF
#define N_CONF 9
#endif /* N_CONF */
#ifndef N_OPT
#define N_OPT 8
#endif /* N_OPT */

double * RHO0 = NULL;
double * U0 = NULL;
double * V0 = NULL;
double * P0 = NULL;
double * X0 = NULL;
double * Y0 = NULL;

int main(int argc, char *argv[])
{
  printf("*********************************************************\n");

  int stat_mkdir = 0, len;
  char add_mkdir[100] = "../SOLUTION/\0";
  DIR * dir_test = NULL;
  strcat(add_mkdir, argv[2]);
  strcat(add_mkdir, "\0");
  dir_test = opendir(add_mkdir);
  if(dir_test != NULL)
    printf("\nOutput directory \"%s\"\n already exists.\n\n", add_mkdir);
  else
  {
    stat_mkdir = mkdir(add_mkdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(stat_mkdir)
    {
      printf("\nOutput directory \"%s\"\n construction failed.\n\n", add_mkdir);
      exit(9);
    }
    else
      printf("\nOutput directory \"%s\"\n constructed.\n\n", add_mkdir);
  }
  closedir(dir_test);
  strcat(add_mkdir, "/rho\0");
  mkdir(add_mkdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  len = strlen(add_mkdir);
  add_mkdir[len-3] = 'm';
  add_mkdir[len-2] = 'e';
  add_mkdir[len-1] = 's';
  add_mkdir[len]   = 'h';
  add_mkdir[len+1] = '\0';
  mkdir(add_mkdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


  int len_prob = 0, i = 0, l = 0, j = 0, k = 0;
  // len_prob is the number of characters in argv[1] EXCLUDING '\0'
  while(argv[1][len_prob] != '\0')
    ++len_prob;
  char addRHO[100];
  char addU[100];
  char addV[100];
  char addP[100];
  char addX[100];
  char addY[100];
  char addCONF[100];
  char addOPT[100];
  char extend[5] = ".txt\0";
  char add[100] = "../DATA/";
  for(i = 0; i < len_prob; ++i)
    add[8+i] = argv[1][i];
  add[8+len_prob] = '/';
  for(i = 0; i < len_prob+9; ++i)
  {
    addRHO[i] = add[i];
    addU[i] = add[i];
    addV[i] = add[i];
    addP[i] = add[i];
    addX[i] = add[i];
    addY[i] = add[i];
    addCONF[i] = add[i];
    addOPT[i] = add[i];
  }
  addRHO[len_prob+9] = 'R';
  addRHO[len_prob+10] = 'H';
  addRHO[len_prob+11] = 'O';
  addU[len_prob+9] = 'U';
  addV[len_prob+9] = 'V';
  addP[len_prob+9] = 'P';
  addX[len_prob+9] = 'X';
  addY[len_prob+9] = 'Y';
  addCONF[len_prob+9] = 'C';
  addCONF[len_prob+10] = 'O';
  addCONF[len_prob+11] = 'N';
  addCONF[len_prob+12] = 'F';
  addOPT[len_prob+9] = 'O';
  addOPT[len_prob+10] = 'P';
  addOPT[len_prob+11] = 'T';
  for(i = 0; i < len_prob; ++i)
  {
    addRHO[len_prob+12+i] = argv[1][i];
    addU[len_prob+10+i] = argv[1][i];
    addV[len_prob+10+i] = argv[1][i];
    addP[len_prob+10+i] = argv[1][i];
    addX[len_prob+10+i] = argv[1][i];
    addY[len_prob+10+i] = argv[1][i];
    addCONF[len_prob+13+i] = argv[1][i];
    addOPT[len_prob+12+i] = argv[1][i];
  }
  for(i = 0; i < 5; ++i)
  {
    addRHO[len_prob+12+len_prob+i] = extend[i];
    addU[len_prob+10+len_prob+i] = extend[i];
    addV[len_prob+10+len_prob+i] = extend[i];
    addP[len_prob+10+len_prob+i] = extend[i];
    addX[len_prob+10+len_prob+i] = extend[i];
    addY[len_prob+10+len_prob+i] = extend[i];
    addCONF[len_prob+13+len_prob+i] = extend[i];
    addOPT[len_prob+12+len_prob+i] = extend[i];
  }
  i = 0;

  double CONFIG[N_CONF];  /* config[0] is the constant of the perfect gas
                           * config[1] is the CFL number
			   * config[2] is the largest value can be seen as zero
			   * config[3] is the first limiter of the slope
			   * config[4] is the second limiter of the slope
			   * config[5] is the first parameter of the monitor function
			   * config[6] is the second parameter of the monitor function
			   * config[7] is the modifier of the mesh redistribution
			   * config[8] is the tolerance of the mesh redistribution
			   */
  double OPT[N_OPT];  /* OPT[0] is the maximal step to compute.
		       * OPT[1] is the time to stop the computation
		       * OPT[2] is the switch of whether keep the inter-data during the computation
		       * OPT[3] is the switch of whether use an adaptive mesh
		       * OPT[4] denote the kind of boundary condition
		       * OPT[5] indicates whether the initial data are the primitive variables [1],
		       *        or the conservative ones [0]
		       * OPT[6] indicates whether we use the smooth derivatives [0],
		       *        or the WENO-type ones in the reconstruction
		       * OPT[7] is the switch of whether use the limiter in the reconstruction
		       */
  int state;
  state = configurate(CONFIG, argv[0], argv[1], addCONF);
  if(state){
    mem_release();
    exit(state);}
  state = optionize(OPT, argv[0], argv[1], addOPT);
  if(state){
    mem_release();
    exit(state);}
  double gamma = CONFIG[0];
  int adp        = (int)OPT[3];
  int inter_data = (int)OPT[2];
  int MaxStp     = (int)OPT[0];
  int vM = inter_data * MaxStp;  // vM == MaxStp or vN = 0
  double scaling = 1.0;

  state = initialize(argv[0], argv[1], addRHO, addU, addP, addX, adp);  /* Firstly we read the initial
									 * data file. The function 
									 * initialize return a point
									 * pointing to the position
									 * of a block of memory
									 * consisting (m+1) variables
									 * of type double.
									 * The value of first of these
									 * variables is m. The
									 * following m variables
									 * are the initial value.
									 */
  if(state){
    mem_release();
    exit(state);}
	
  int m = (int)RHO0[0];  /* m is the number of initial value
			  * as well as the number of grids.
			  * As m is frequently use to
			  * represent the number of grids,
			  * we do not use the name such as
			  * num_grid here to correspond to
			  * notation in the math theory.
			  */


  runList runhist;
  init_runList(&runhist);

  
  double * rho[vM+1];
  for(k = 0; k <= vM; ++k)
  {
    rho[k] = (double *)malloc(m * sizeof(double));
    if(rho[k] == NULL)
    {
      for(i = 0; i < k; ++i)
      {
	free(rho[i]);
	rho[i] = NULL;
      }
      free(RHO0);
      free(U0);
      free(P0);
      RHO0 = NULL;
      U0 = NULL;
      P0 = NULL;
      printf("NOT enough memory! RHO[%d]\n", k);
      exit(5);
    }
  }
  double * u[vM+1];
  for(k = 0; k <= vM; ++k)
  {
    u[k] = (double *)malloc(m * sizeof(double));
    if(u[k] == NULL)
    {
      for(i = 0; i < k; ++i)
      {
	free(u[i]);
	u[i] = NULL;
      }
      for(i = 0; i <= vM; ++i)
      {
	free(rho[i]);
	rho[i] = NULL;
      }
      free(RHO0);
      free(U0);
      free(P0);
      RHO0 = NULL;
      U0 = NULL;
      P0 = NULL;
      printf("NOT enough memory! U[%d]\n", k);
      exit(5);
    }
  }



  double * p[vM+1];
  p[0] = P0 + 1;
  for(k = 0; k <= vM; ++k)
  {
    p[k] = (double *)malloc(m * sizeof(double));
    if(p[k] == NULL)
    {
      for(i = 0; i < k; ++i)
      {
	free(p[i]);
	p[i] = NULL;
      }
      for(i = 0; i <= vM; ++i)
      {
	free(rho[i]);
	rho[i] = NULL;
	free(u[i]);
	u[i] = NULL;
      }
      free(RHO0);
      free(U0);
      free(P0);
      RHO0 = NULL;
      U0 = NULL;
      P0 = NULL;
      printf("NOT enough memory! P[%d]\n", k);
      exit(5);
    }
  }


  for(j = 0; j < m; ++j)
  {
    rho[0][j] = RHO0[j+1];
    u[0][j] =   U0[j+1];
    p[0][j] =   P0[j+1];
  }


  double h;
  // vm = m+1;
  int vm = adp*m+1;
  double *x[vM+1], *xc[vM+1];


  for(k = 0; k <= vM; ++k)
  {
    x[k] = (double *)malloc(sizeof(double)*vm);
    if(x[k] == NULL)
    {
      for(i = 0; i <= vM; ++i)
      {
	free(rho[i]);
	rho[i] = NULL;
	free(u[i]);
	u[i] = NULL;
	free(p[i]);
	p[i] = NULL;
      }
      for(i = 0; i < k; ++i)
      {
	free(x[i]);
	x[i] = NULL;
      }
      free(RHO0);
      free(U0);
      free(P0);
      free(X0);
      RHO0 = NULL;
      U0 = NULL;
      P0 = NULL;
      X0 = NULL;
      printf("NOT enough memory! P[%d]\n", k);
      exit(5);
    }
  }
  for(k = 0; k <= vM; ++k)
  {
    xc[k] = (double *)malloc(sizeof(double)*m);
    if(xc[k] == NULL)
    {
      for(i = 0; i <= vM; ++i)
      {
	free(rho[i]);
	rho[i] = NULL;
	free(u[i]);
	u[i] = NULL;
	free(p[i]);
	p[i] = NULL;
	free(x[i]);
	x[i] = NULL;
      }
      for(i = 0; i < k; ++i)
      {
	free(xc[i]);
	xc[i] = NULL;
      }
      free(RHO0);
      free(U0);
      free(P0);
      free(X0);
      RHO0 = NULL;
      U0 = NULL;
      P0 = NULL;
      X0 = NULL;
      printf("NOT enough memory! P[%d]\n", k);
      exit(5);
    }
  }

  if(adp)
    for(j = 0; j < vm; ++j)
      x[0][j] = X0[j+1]*scaling;
  else
    for(k = 0; k <= vM; ++k)
    {
      free(xc[k]);
      xc[k] = NULL;
    }
  h = X0[1]*scaling;
  OPT[1] = OPT[1]*scaling;


  free(RHO0);
  free(U0);
  free(P0);
  free(X0);
  RHO0 = NULL;
  U0 = NULL;
  P0 = NULL;
  X0 = NULL;


  int K = 0;
  char scheme[100];
  //*
  K = GRP_HWENO_fix(CONFIG, OPT, m, h, rho, u, p, &runhist, scheme);
  file_write_trouble(m, K, &runhist, argv[2]);/*/
  K = GRP_fix(CONFIG, OPT, m, h, rho, u, p, &runhist, scheme);//*/


  int vvM = vM;
  if(vvM > K)
    vvM = K;
  int start = vvM-5;
  if(start < 0)
    start = 0;
  start = 0;

  /*
   * vM == MaxStp or vM = 0
   * K no greater than MaxStp
   *
   * if adp
   *   vvM = K <= MaxStp
   * if fix
   *   vvM = 0
   *
   * runhist[1 -- K]
   * data[0 -- vvM]
   */

  file_write_log(m, 1, K, scaling, CONFIG, OPT, &runhist, scheme, argv[1], argv[2]);
  file_write_data(m, start, vvM, rho, "rho", argv[2]);
  file_write_data(m, start, vvM,   u, "u__", argv[2]);
  file_write_data(m, start, vvM,   p, "p__", argv[2]);
  if(adp)
  {
    file_write_data(m, start, vvM,  xc, "xc_", argv[2]);
    file_write_data(m+1, start, vvM, x, "x__", argv[2]);
  }



  delete_runList(&runhist);
  for(k = 0; k <= vM; ++k)
  {
    free(rho[k]);
    free(u[k]);
    free(p[k]);
    rho[k] = NULL;
    u[k] = NULL;
    p[k] = NULL;
  }
  if(adp)
  {
    for(k = 0; k <= vM; ++k)
    {
      free(xc[k]);
      xc[k] = NULL;
    }
    for(k = 0; k <= vM; ++k)
    {
      free(x[k]);
    }
  }



  printf("\n");
  return 0;
}
