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

#include "file_io.h"
#ifndef L_STR
#include "file_io_local.h"
#endif
#include "solver.h"
#include "reconstruction.h"






int main(int argc, char *argv[])
{
  int i = 0, l = 0, j = 0, k = 0, it;
  char err_msg[L_STR];
  int read_state, err_code = 0;
  double CONFIG[N_CONF];
  int already_read[N_CONF];
  char ITEM[N_CONF][L_STR];


/* CONFIG[0]  is the constant of the perfect gas
 * CONFIG[1]  is the CFL number
 * CONFIG[2]  is the largest value can be seen as zero
 * CONFIG[3]  is the first limiter of the slope
 * CONFIG[4]  is the second limiter of the slope
 * CONFIG[5]  is the first parameter of the monitor function
 * CONFIG[6]  is the second parameter of the monitor function
 * CONFIG[7]  is the modifier of the mesh redistribution
 * CONFIG[8]  is the tolerance of the mesh redistribution
 * CONFIG[9]  is the parameter of identifying the discontinuities
 * CONFIG[10] is the parameter of the thickness in the THINC reconstruction
 * CONFIG[11] is the maximal step to compute.
 * CONFIG[12] is the switch of whether use an adaptive mesh
 * CONFIG[13] denote the kind of boundary condition
 * CONFIG[14] indicates whether the initial data are the primitive
 *            variables [1], or the conservative ones [0]
 * CONFIG[15] indicates whether we use the smooth derivatives [0],
 *            or the WENO-type ones in the reconstruction
 * CONFIG[16] is the switch of whether use the limiter in the reconstruction
 * CONFIG[17] is the switch of whether use the characteristic decomposition
 * CONFIG[18] is the scaling
 */
  read_state = configurate(CONFIG, ITEM, already_read, err_msg, "CONF\0", argv[1]);
  if(read_state)
  {
    printf("%s", err_msg);
    err_code = read_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, read_state-(err_code<<bit_shift));
    return read_state;
  }
  display_config(CONFIG, argv[1]);
  int adp        = (int)CONFIG[12];
  double scaling = CONFIG[18];


  int nInitValue = 5;
  char addInitValue[nInitValue][L_STR];
  strcpy(addInitValue[0], "RHO\0");
  strcpy(addInitValue[1], "U\0\0");
  strcpy(addInitValue[2], "P\0");
  strcpy(addInitValue[3], "X\0");
  strcpy(addInitValue[4], "TIME\0");
  int sizeInitValue[nInitValue];
  realArray InitValue[nInitValue];
  for(it = 0; it < nInitValue; ++it)
    InitValue[it].n_box = 0;


  read_state = initialize(nInitValue, InitValue, sizeInitValue, err_msg, addInitValue, argv[1]);
  if(read_state)
  {
    printf("%s", err_msg);
    err_code = read_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, read_state-(err_code<<bit_shift));
    return read_state;
  }
  int m = sizeInitValue[0];




  
  double * rho;
  rho = (double *)malloc(m * sizeof(double));
  if(rho == NULL)
  {
    for(k = 0; k < nInitValue; ++k)
      delete_realArray(InitValue +k);
    printf("NOT enough memory! RHO\n");
    exit(11);
  }
  double * u;
  u = (double *)malloc(m * sizeof(double));
  if(u == NULL)
  {
    free(rho);
    rho = NULL;
    for(k = 0; k < nInitValue; ++k)
      delete_realArray(InitValue +k);
    printf("NOT enough memory! U\n");
    exit(12);
  }
  double * p;
  p = (double *)malloc(m * sizeof(double));
  if(p == NULL)
  {
    free(rho);
    rho = NULL;
    free(u);
    u = NULL;
    for(k = 0; k < nInitValue; ++k)
      delete_realArray(InitValue +k);
    printf("NOT enough memory! P\n");
    exit(13);
  }

  int vm = adp*m+1;
  double * x, * xc;
  double h;

  x = (double *)malloc(sizeof(double)*vm);
  if(x == NULL)
  {
    free(rho);
    rho = NULL;
    free(u);
    u = NULL;
    free(p);
    p = NULL;
    for(k = 0; k < nInitValue; ++k)
      delete_realArray(InitValue +k);
    printf("NOT enough memory! X\n");
    exit(16);
  }
  xc = (double *)malloc(sizeof(double)*m);
  if(xc == NULL)
  {
    free(rho);
    rho = NULL;
    free(u);
    u = NULL;
    free(p);
    p = NULL;
    free(x);
    x = NULL;
    for(k = 0; k < nInitValue; ++k)
      delete_realArray(InitValue +k);
    printf("NOT enough memory! XC\n");
    exit(17);
  }
  double *TIME;
  int nTIME = sizeInitValue[4];
  TIME = (double *)malloc(sizeof(double)*nTIME);



  double * DATA_OUT[nInitValue];
  DATA_OUT[0] = rho;
  DATA_OUT[1] = u;
  DATA_OUT[2] = p;
  DATA_OUT[3] = x;
  DATA_OUT[4] = TIME;
  
  for(it = 0; it < nInitValue; ++it)
    for(j = 0; j < sizeInitValue[it]; ++j)
      find_cargo_realArray(DATA_OUT[it]+j, j, InitValue+it);
  
  if(adp)
    for(j = 0; j < vm; ++j)
      x[j] = x[j];
  else
  {
    free(xc);
    xc = NULL;
    h = x[0];
  }


  for(k = 0; k < nInitValue; ++k)
    delete_realArray(InitValue +k);



  int len;
  for(it = 0; it < nInitValue; ++it)
  {
    len = 0;
    while(addInitValue[it][len] != '\0')
      ++len;
    for(i = 0; i < len; ++i)
      addInitValue[it][i] += 32;
    addInitValue[it][len] = '\0';
  }


  int output_flag[nInitValue], output_idx[nInitValue], output_state;
  output_flag[0] = 1; //rho
  output_flag[1] = 1; //u
  output_flag[2] = 1; //p
  output_flag[3] = 0; //x
  output_flag[4] = 0; //time
  for(it = 0; it < nInitValue; ++it)
    output_idx[it] = nTIME-1;

  option OPT;
  OPT.gamma     = CONFIG[0];  // the constant of the perfect gas
  OPT.CFL       = CONFIG[1];  // CFL number
  OPT.eps       = CONFIG[2];  // the largest value could be treat as zero
  OPT.tol       = CONFIG[8];
  OPT.MaxStp    = (int)(CONFIG[11]);  // the number of time steps
  OPT.TIME      = DATA_OUT[4];
  OPT.nTIME     = sizeInitValue[4];

  OPT.bod       = (int)CONFIG[13];
  OPT.Primative = (int)CONFIG[14];
  OPT.Deri      = (int)CONFIG[15];
  OPT.Limiter   = (int)CONFIG[16];
  OPT.Decomp    = (int)CONFIG[17];

  OPT.alp1      = CONFIG[3];
  OPT.alp2      = CONFIG[4];
  OPT.bet1      = CONFIG[5];
  OPT.bet2      = CONFIG[6];
  OPT.modifier  = CONFIG[7];
  OPT.threshold = CONFIG[9];
  OPT.thickness = CONFIG[10];

  OPT.output_flag = output_flag;
  OPT.nInitValue = nInitValue;
  OPT.sizeInitValue = sizeInitValue;
  OPT.addInitValue = (char **) malloc(sizeof(char *) * nInitValue);
  for(it = 0; it < nInitValue; ++it)
    OPT.addInitValue[it] = addInitValue[it];




  runHist runhist;
  init_runHist(&runhist);
  int K = 0;
  char scheme[L_STR];
  char version[L_STR] = "dev\0";
  char add_mkdir[L_STR+L_STR];
  strcpy(add_mkdir, version);
  printf("The present version is [%s]\n", version);
  //K = GRP2_fix(OPT, m, h, rho, u, p, &runhist, add_mkdir, argv[2]);
  K = GRP4_HWENO5_fix(OPT, m, h, rho, u, p, &runhist, add_mkdir, argv[2]);
  //K = RF4_WENO5_fix(OPT, m, h, rho, u, p, &runhist, add_mkdir, argv[2]);

  
  output_state = DATA_OUTPUT(err_msg, nInitValue, addInitValue, sizeInitValue, DATA_OUT, output_idx, output_flag, add_mkdir);
  if(output_state)
  {
    printf("%s", err_msg);
    err_code = output_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, output_state-(err_code<<bit_shift));
    delete_runHist(&runhist);
    free(rho);
    free(u);
    free(p);
    free(x);
    rho = NULL;
    u = NULL;
    p = NULL;
    x = NULL;
    if(adp)
    {
      free(xc);
      xc = NULL;
    }
    free(TIME);
    return output_state;
  }
  printf("DATA OUT.\n");

  
  output_state = LOG_OUTPUT(err_msg, &runhist, ITEM, already_read, CONFIG, m, 1, K, scheme, version, argv[1], add_mkdir);
  printf("LOG OUT.\n");


  delete_runHist(&runhist);
  
  printf("RUN DELETE.\n");
  free(rho);
  free(u);
  free(p);
  free(x);
  rho = NULL;
  u = NULL;
  p = NULL;
  x = NULL;
  if(adp)
  {
    free(xc);
    xc = NULL;
  }
  free(TIME);
  printf("\n");
  return 0;
}
