#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "file_io.h"



#ifndef N_CONF
#define N_CONF 9
#endif /* N_CONF */
#ifndef N_OPT
#define N_OPT 8
#endif /* N_OPT */


extern double * RHO0;
extern double * U0;
extern double * V0;
extern double * P0;
extern double * X0;
extern double * Y0;


void mem_release()
{
  if(RHO0)
    free(RHO0);
  if(U0)
    free(U0);
  if(P0)
    free(P0);
  if(X0)
    free(X0);
}



int vec_read(char * add, int COLUMN, int err_code)
{
  int line, column, num, file_read_state;
  FILE * fp_data;
  double * data;
  char source;

  if(err_code == 10)
    source = 'R';
  else if(err_code == 20)
    source = 'U';
  else if(err_code == 40)
    source = 'P';
  else if(err_code == 50)
    source = 'X';
  else{
    printf("err_code = %d", err_code);
    return 100;}

  if((fp_data = fopen(add, "r")) == 0){
    printf("Cannot open initial data file %s!\n", add);
    return(err_code);}

  line = data_pre_read_line(fp_data, err_code);
  if(line < 0){
    fclose(fp_data);
    return line;}
  else if(line > 1){
    printf("The initial data [%c] has %d lines!\n", source, line);
    fclose(fp_data);
    return(err_code + 2);}
  fseek(fp_data, 0L, SEEK_SET);

  column = data_pre_read_column(fp_data, err_code);
  if(column < 0){
    fclose(fp_data);
    return column;}
  if((column != COLUMN) && (err_code != 10)){
    printf("[%c] UNEQUAL! column_u=%d\tcolumn_rho=%d\n", source, column, COLUMN);
    fclose(fp_data);
    return(err_code + 3);}
  fseek(fp_data, 0L, SEEK_SET);





  if(err_code == 10){
    RHO0 = (double *)malloc((column + 1) * sizeof(double));
    data = RHO0;}
  else if(err_code == 20){
    U0 = (double *)malloc((column + 1) * sizeof(double));
    data = U0;}
  else if(err_code == 30){
    V0 = (double *)malloc((column + 1) * sizeof(double));
    data = V0;}
  else if(err_code == 40){
    P0 = (double *)malloc((column + 1) * sizeof(double));
    data = P0;}
  else if(err_code == 50){
    X0 = (double *)malloc((column + 1) * sizeof(double));
    data = X0;}
  else{
    Y0 = (double *)malloc((column + 1) * sizeof(double));
    data = Y0;}

  data[0] = (double)column;
  file_read_state = data_read(fp_data, data+1, num);
  fclose(fp_data);
  data = NULL;

  if(file_read_state)
  {
    if(file_read_state == column)
      printf("Error during file reading! [%s]\n", add);
    else
      printf("\nThe %dth entry in the file [%s] is not a number.\n", file_read_state, add);
    return(err_code + 5);
  }

  return 0;
}



/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
int initialize(char * scheme, char * label, char * addrho, char * addu, char * addp, char * addx, int adp)
{
  int read_state;
  int COLUMN;
  int file_read_state;


  read_state = vec_read(addrho, 0, 10);
  if(read_state){
    return read_state;}
  COLUMN = RHO0[0];

  read_state = vec_read(addu, COLUMN, 20);
  if(read_state){
    return read_state;}

  read_state = vec_read(addp, COLUMN, 40);
  if(read_state){
    return read_state;}

  if(adp)
  {
    read_state = vec_read(addx, COLUMN+1, 50);
    if(read_state){
      return read_state;}
  }
  else
  {
    read_state = vec_read(addx, 1, 50);
    if(read_state){
      return read_state;}
  }
  printf("[%s | %s] initialized, column=%d.\n\n\n", scheme+2, label, COLUMN);

  return 0;
}


/* This function read the configuration data file,
 * and store the configuration data in the array
 * "CONFIG".
 * CONFIG[0] is the constant of perfect gas
 * CONFIG[1] is the CFL number
 * CONFIG[2] is the largest value can be seen as zero
 * CONFIG[3] is the first limiter of the slope
 * CONFIG[4] is the second limiter of the slope
 * CONFIG[5] is the first parameter of the monitor function
 * CONFIG[6] is the second parameter of the monitor function
 * CONFIG[7] is the modifier of the mesh redistribution
 * CONFIG[8] is the tolerance of the mesh redistribution
 */
int configurate(double * CONFIG, char * scheme, char * label, char * add)
{
  FILE * fp_data;
  int n_conf, state, k, j, compare[N_CONF], sgn = 1;
  char * str_conf[N_CONF*2];
  int const item_L = 101;

  char ITEM[N_CONF][10] = {"gamma", "CFL", "eps", "alp1", "alp2", "bet1", "bet2", "modifier", "tol"};
  double up_bound[2][N_CONF] = {{0.0, 1.0, 0.1, 2.0, 2.0, 0.0, 0.0, 1.0, 0.1},
				{0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}};
  double low_bound[N_CONF] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


  for(k = 0; k < N_CONF*2; ++k)
    str_conf[k] = malloc(item_L * sizeof(char));
  for(k = 0; k < N_CONF; ++k)
    compare[k] = 0;

  //open the configuration data file
  //printf("%s will open the cinfiguration data file: ", name);
  //printf("%s\n\n", add);
  if((fp_data = fopen(add, "r")) == 0)
  {
    printf("Cannot open configuration data file: %s!\n", add);
    return 70;
  }

  n_conf = string_pre_read_line(fp_data, 70);
  if(n_conf < 0){
    fclose(fp_data);
    return n_conf;}
  if(n_conf != N_CONF){
    printf("%d, %d", n_conf, N_CONF);
    fclose(fp_data);
    return 72;}
  fseek(fp_data, 0L, SEEK_SET);
  n_conf = string_pre_read_column(fp_data, 70);
  if(n_conf < 0){
    fclose(fp_data);
    return n_conf;}
  if(n_conf != 2){
    printf("%d, 2", n_conf);
    fclose(fp_data);
    return 73;}
  fseek(fp_data, 0L, SEEK_SET);

  state =  string_read(fp_data, str_conf, N_CONF*2, item_L-1);
  fclose(fp_data);
  if(state)
  {
    if(state > 0)
      printf("Too many items in the configuration file! [%d | %d]\n", state, N_CONF*2);
    else
      printf("Too many charactors in the %d-th item in the configuration data file\n", -state);
    return(75);
  }



  for(k = 0; k < N_CONF; ++k)
  {
    for(j = 0; j < N_CONF; ++j)
      if(!strcmp(str_conf[2*k], ITEM[j]))
	break;

    if(j == N_CONF)
    {
      printf("Illigal item [%d: %s] in the configuration file.\n", 2*k, str_conf[2*k]);
      return 76;
    }
    if(compare[j])
    {
      printf("Doubled item [%s] in the configuration file.\n", str_conf[2*k]);
      return 76;
    }

    if(!(sgn = format_string(str_conf[2*k+1])))
    {
      printf("The %d-th item in the configuration file is not a number.\n", k+1);
      return 77;
    }
    CONFIG[j] = sgn*str2num(str_conf[2*k+1]);
    compare[j] = 1;

    if((up_bound[1][j]*CONFIG[j] > up_bound[0][j]) || (CONFIG[j] < low_bound[j]))
    {
      printf("Illigal value [%g] of the item [%s].\n", CONFIG[j], str_conf[2*k]);
      return 78;
    }
  }

  for(k = 0; k < N_CONF*2; ++k)
    free(str_conf[k]);


  printf("[%s | %s] configurated:\n", scheme+2, label);
  printf("  gamma  = %g\n", CONFIG[0]);
  printf("  CFL    = %g\n", CONFIG[1]);
  printf("  eps    = %g\n", CONFIG[2]);
  printf("  alpha1 = %g\tthe first  limiter of the slope\n", CONFIG[3]);
  printf("  alpha2 = %g\tthe second limiter of the slope\n", CONFIG[4]);
  printf("  beta1  = %g\tthe first  parameter of the monitor function\n", CONFIG[5]);
  printf("  beta2  = %g\tthe second parameter of the monitor function\n", CONFIG[6]);
  printf("  modfr  = %g\n", CONFIG[7]);
  printf("  tol    = %g\n", CONFIG[8]);

  return 0;
}

/* OPT[0] is the maximal step to compute.
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
int optionize(double * OPT, char * scheme, char * label, char * add)
{
  FILE * fp_data;
  int n_opt, state, k, j, compare[N_OPT], sgn = 1;
  char * str_opt[N_OPT*2];
  int const item_L = 101;

  char ITEM[N_OPT][10] = {"MaxStp", "Time", "InterData", "Adp", "boundary", "Riemann", "WENOD", "limiter"};


  for(k = 0; k < N_OPT*2; ++k)
    str_opt[k] = malloc(item_L * sizeof(char));
  for(k = 0; k < N_OPT; ++k)
    compare[k] = 0;

  if((fp_data = fopen(add, "r")) == 0)
  {
    printf("Cannot open option data file: %s!\n", add);
    return 80;
  }

  n_opt = string_pre_read_line(fp_data, 80);
  if(n_opt < 0){
    fclose(fp_data);
    return n_opt;}
  if(n_opt != N_OPT){
    printf("n_opt=%d, N_OPT=%d", n_opt, N_OPT);
    fclose(fp_data);
    return 82;}
  fseek(fp_data, 0L, SEEK_SET);
  n_opt = string_pre_read_column(fp_data, 80);
  if(n_opt < 0){
    fclose(fp_data);
    return n_opt;}
  if(n_opt != 2){
    printf("%d, 2", n_opt);
    fclose(fp_data);
    return 83;}
  fseek(fp_data, 0L, SEEK_SET);

  state =  string_read(fp_data, str_opt, N_OPT*2, item_L-1);
  fclose(fp_data);
  if(state)
  {
    if(state > 0)
      printf("Too many items in the option file! [%d | %d]\n", state, N_OPT*2);
    else
      printf("Too many charactors in the %d-th item in the option data file\n", -state);
    return(85);
  }


  for(k = 0; k < N_OPT; ++k)
  {
    for(j = 0; j < N_OPT; ++j)
      if(!strcmp(str_opt[2*k], ITEM[j]))
	break;

    if(j == N_OPT)
    {
      printf("Illigal item [%d: %s] in the option file.\n", 2*k, str_opt[2*k]);
      return 86;
    }
    if(compare[j])
    {
      printf("Doubled item [%s] in the option file.\n", str_opt[2*k]);
      return 86;
    }

    if(!(sgn = format_string(str_opt[2*k+1])))
    {
      printf("The %d-th item in the option file is not a number.\n", k+1);
      return 87;
    }
    OPT[j] = sgn*str2num(str_opt[2*k+1]);
    compare[j] = 1;
  }

  for(k = 0; k < N_OPT*2; ++k)
    free(str_opt[k]);


  if(OPT[0] < 0.0){
    printf("The maximum step should be positive. [%g]\n", OPT[0]);
    return 88;}
  printf("  MaxStp = %d\n", (int)OPT[0]);
  if(OPT[1] < 0.0){
    printf("The TIME should be positive. [%g]\n", OPT[1]);
    return 88;}
  printf("  TIME   = %g\n", OPT[1]);
  //printf("\n");
  int inter_data = (int)OPT[2];
  if(inter_data)
    printf("  inter-data     kept\n");
  else
    printf("  inter-data NOT kept\n");

  int adp = (int)OPT[3];
  if(adp)
    printf("  ADAPTIVE       mesh\n");
  else
    printf("  FIXED          mesh\n");

  int bod = (int)OPT[4];
  if(bod < 0)
    printf("  REFLECTION     boundary condition\n");
  else if(bod > 0)
    printf("  INFINITY       boundary condition\n");
  else
    printf("  PERIODIC       boundary condition\n");

  int Riemann = (int)OPT[5];
  if(Riemann)
    printf("  PRIMITIVE      variables\n");
  else
    printf("  CONSERVATION   variables\n");

  int WENOD = (int)OPT[6];
  if(WENOD)
    printf("  WENO-type      direvatives\n");
  else
    printf("  SMOOTH         direvatives\n");

  int limiter = (int)OPT[7];
  if(limiter)
    printf("  limiter        used\n");
  else
    printf("  limiter    NOT used\n");
  printf("\n");

  return 0;
}
