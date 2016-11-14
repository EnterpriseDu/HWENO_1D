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
 * CONFIG[12] is the time to stop the computation
 * CONFIG[13] is the switch of whether use an adaptive mesh
 * CONFIG[14] denote the kind of boundary condition
 * CONFIG[15] indicates whether the initial data are the primitive
 *            variables [1], or the conservative ones [0]
 * CONFIG[16] indicates whether we use the smooth derivatives [0],
 *            or the WENO-type ones in the reconstruction
 * CONFIG[17] is the switch of whether use the limiter in the reconstruction
 * CONFIG[18] is the switch of whether use the characteristic decomposition
 * CONFIG[19] is the scaling
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
  int adp        = (int)CONFIG[13];
  double scaling = CONFIG[19];


  int nInitValue = 4;
  char addInitValue[nInitValue][L_STR];
  strcpy(addInitValue[0], "RHO\0");
  strcpy(addInitValue[1], "U\0\0");
  strcpy(addInitValue[2], "P\0");
  strcpy(addInitValue[3], "X\0");
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



  double * DATA_OUT[nInitValue];
  DATA_OUT[0] = rho;
  DATA_OUT[1] = u;
  DATA_OUT[2] = p;
  DATA_OUT[3] = x;
  
  for(it = 0; it < nInitValue; ++it)
    for(j = 0; j < sizeInitValue[it]; ++j)
      find_cargo_realArray(DATA_OUT[it]+j, j, InitValue+it);
  
  if(adp)
    for(j = 0; j < vm; ++j)
      x[j] = x[j] * scaling;
  else
  {
    free(xc);
    xc = NULL;
    h = x[0]*scaling;
  }
  CONFIG[12] = CONFIG[12]*scaling;


  for(k = 0; k < nInitValue; ++k)
    delete_realArray(InitValue +k);






  runHist runhist;
  init_runHist(&runhist);
  int K = 0;
  char scheme[L_STR];
  char version[L_STR] = "dev\0";
  char add_mkdir[L_STR+L_STR];
  strcpy(add_mkdir, version);
  printf("The present version is [%s]\n", version);
  K = Exct2_2_fix(CONFIG, m, h, rho, u, p, &runhist, add_mkdir, argv[2]);
  //K = GRP2_fix(CONFIG, m, h, rho, u, p, &runhist, add_mkdir, argv[2]);


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
  for(it = 0; it < nInitValue; ++it)
    output_idx[it] = 0;
  
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
  printf("\n");
  return 0;
}
