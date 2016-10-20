#include <stdio.h>
#include <stdlib.h>


#include "file_io.h"
#ifndef L_STR
#include "file_io_local.h"
#endif



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
int configurate(double CONFIG[N_CONF], char ITEM[N_CONF][L_STR], int already_read[N_CONF], char * err_msg, char * addCONF, char * prob)
{
  int err_code = 01;
  char extend[5] = ".txt\0";
  char add[L_STR+L_STR];
  int len_prob = 0, len_CONF = 0, i, j, k, l, count, idx;
  FILE * fp_data;
  strcpy(add, "../DATA/\0");

  while(prob[len_prob] != '\0')
    ++len_prob;
  while(addCONF[len_CONF] != '\0')
    ++len_CONF;

  idx = 8;
  for(i = 0; i < len_prob; ++i)
    add[idx++] = prob[i];
  add[idx++] = '/';
  for(i = 0; i < len_CONF; ++i)
    add[idx++] = addCONF[i];
  for(i = 0; i < len_prob; ++i)
    add[idx++] = prob[i];
  for(i = 0; i < 5; ++i)
    add[idx++] = extend[i];

  for(i = 0; i < N_CONF; ++i)
    already_read[i] = 0;
  //the value of already_read[idx] is the position of ITEM[idx] occured in the file
  int read_flag = 0;//0 for an item, 1 for a value
  double up_bound[2][N_CONF];
  double low_bound[N_CONF];
  double default_value[N_CONF];
  strcpy(ITEM[0], "gamma\0");
  up_bound[0][0]   = 0.0;
  up_bound[1][0]   = 0.0;
  low_bound[0]     = 1.0;
  default_value[0] = 1.4;
  strcpy(ITEM[1], "CFL\0");
  up_bound[0][1]   = 1.0;
  up_bound[1][1]   = 1.0;
  low_bound[1]     = 0.0;
  default_value[1] = 0.45;
  strcpy(ITEM[2], "eps\0");
  up_bound[0][2]   = 1e-5;
  up_bound[1][2]   = 1.0;
  low_bound[2]     = 0.0;
  default_value[2] = 1e-9;
  strcpy(ITEM[3], "alp1\0");
  up_bound[0][3]   = 1.0;
  up_bound[1][3]   = 1.0;
  low_bound[3]     = 0.0;
  default_value[3] = 1.8;
  strcpy(ITEM[4], "alp2\0");
  up_bound[0][4]   = 2.0;
  up_bound[1][4]   = 1.0;
  low_bound[4]     = 0.0;
  default_value[4] = 0.9;
  strcpy(ITEM[5], "bet1\0");
  up_bound[0][5]   = 0.0;
  up_bound[1][5]   = 0.0;
  low_bound[5]     = 0.0;
  default_value[5] = 20.0;
  strcpy(ITEM[6], "bet2\0");
  up_bound[0][6]   = 0.0;
  up_bound[1][6]   = 0.0;
  low_bound[6]     = 0.0;
  default_value[6] = 100.0;
  strcpy(ITEM[7], "modifier\0");
  up_bound[0][7]   = 1.0;
  up_bound[1][7]   = 1.0;
  low_bound[7]     = 0.0;
  default_value[7] = 0.05;
  strcpy(ITEM[8], "tol\0");
  up_bound[0][8]   = 1e-4;
  up_bound[1][8]   = 1.0;
  low_bound[8]     = 0.0;
  default_value[8] = 1e-6;
  strcpy(ITEM[9], "threshold\0");
  up_bound[0][9]   = 0.0;
  up_bound[1][9]   = 0.0;
  low_bound[9]     = 0.0;
  default_value[9] = 15.0;
  strcpy(ITEM[10], "thickness\0");
  up_bound[0][10]   = 0.0;
  up_bound[1][10]   = 0.0;
  low_bound[10]     = 0.0;
  default_value[10] = 0.8;
  strcpy(ITEM[11], "MaxStp\0");
  up_bound[0][11]   = 0.0;
  up_bound[1][11]   = 0.0;
  low_bound[11]     = 0.0;
  default_value[11] = 0.0;
  strcpy(ITEM[12], "Time\0");
  up_bound[0][12]   = 0.0;
  up_bound[1][12]   = 0.0;
  low_bound[12]     = 0.0;
  default_value[12] = 0.0;
  strcpy(ITEM[13], "Adp\0");
  up_bound[0][13]   = 0.0;
  up_bound[1][13]   = 0.0;
  low_bound[13]     = 0.0;
  default_value[13] = 0.0;
  strcpy(ITEM[14], "Boundary\0");
  up_bound[0][14]   = 0.0;
  up_bound[1][14]   = 0.0;
  low_bound[14]     = 0.0;
  default_value[14] = 1.0;
  strcpy(ITEM[15], "Primative\0");
  up_bound[0][15]   = 0.0;
  up_bound[1][15]   = 0.0;
  low_bound[15]     = 0.0;
  default_value[15] = 0.0;
  strcpy(ITEM[16], "Deri\0");//zero stands for using smooth interpolation
  up_bound[0][16]   = 0.0;
  up_bound[1][16]   = 0.0;
  low_bound[16]     = 0.0;
  default_value[16] = 0.0;
  strcpy(ITEM[17], "Limiter\0");
  up_bound[0][17]   = 0.0;
  up_bound[1][17]   = 0.0;
  low_bound[17]     = 0.0;
  default_value[17] = 0.0;
  strcpy(ITEM[18], "Decomp\0");
  up_bound[0][18]   = 0.0;
  up_bound[1][18]   = 0.0;
  low_bound[18]     = 0.0;
  default_value[18] = 0.0;
  strcpy(ITEM[19], "scaling\0");
  up_bound[0][19]   = 0.0;
  up_bound[1][19]   = 0.0;
  low_bound[19]     = 1.0;
  default_value[19] = 1.0;

  char comment = '#';
  int n_digit_real = 15;
  char digit_real[] = {'0','1','2','3','4','5','6','7','8','9','-','+','.','e','E'};
  int n_space = 4;
  char space[] = {' ', '\t', '\n', '\r'};
  int read_state;
  char buffer[L_BUFF], number[L_BUFF], end_mark;
  int sign;
  Text text;


  if((fp_data = fopen(add, "r")) == 0)
  {
    sprintf(err_msg, "Cannot open configuration data file: %s!\n", add);
    return err_code;
  }

  init_Text(&text);
  read_state = string_read(&text, &end_mark, fp_data, buffer, '#', EOF, space, n_space, err_msg);
  if(read_state)
  {
    //strcat(err_msg, sprintf(" (the %d-th item)\n", read_state));
    sprintf(err_msg, "%s (the %d-th item)\n", err_msg, read_state);
    delete_Text(&text);
    return err_code + ((01) << bit_shift);
  }

  text.current = text.head;
  count = 0;
  read_flag = 0;//0 for an item, 1 for a value
  while(text.current)
  {
    ++count;
    if(!read_flag)
    {
      for(idx = 0; idx < N_CONF; ++idx)
	if(!strcmp(text.current->words, ITEM[idx]))
	  break;
      if(idx == N_CONF)
      {
	sprintf(err_msg, "Expecting a configuration item, while the %d-th entry '%s' is not one.\n", count, text.current->words);
	delete_Text(&text);
	return err_code + ((02) << bit_shift);
      }
      if(already_read[idx])
      {
	sprintf(err_msg, "Item '%s' double defined. %d-th and %d-th.\n", ITEM[idx], already_read[idx], count);
	delete_Text(&text);
	return err_code + ((03) << bit_shift);
      }

      read_flag = 1;
      already_read[idx] = count;
      //the value of already_read[idx] is the position of ITEM[idx] occured in the file
    }
    else
    {
      strcpy(number, text.current->words);
      sign = is_real(number, digit_real, n_digit_real);
      if(!sign)
      {
	sprintf(err_msg, "Expecting a real number for the value of '%s(%d)' but the %d-th entry is not.\n", ITEM[idx], already_read[idx], count);
	delete_Text(&text);
	return err_code + ((04) << bit_shift);
      }

      CONFIG[idx] = sign*str2real(number);

      if((up_bound[1][idx]*CONFIG[idx] > up_bound[0][idx]) || (CONFIG[idx] < low_bound[idx]))
      {
	sprintf(err_msg, "Illigal value [%g] for %d-th the item '%s'.\n", CONFIG[idx], idx, ITEM[idx]);
	delete_Text(&text);
	return err_code + ((05) << bit_shift);
      }

      //printf("%d, %s: %g\n", idx, ITEM[idx], CONFIG[idx]);

      read_flag = 0;
    }

    text.current = text.current->next;
  }


  for(i = 0; i < N_CONF; ++i)
    if(!already_read[i])
    {
      if(!default_value[i])
      {
	sprintf(err_msg, "Found no value for '%s'.\n", ITEM[i]);
	delete_Text(&text);
	return err_code + ((06) << bit_shift);
      }
      else
      {
	CONFIG[i] = default_value[i];
	already_read[i] = -2;
      }
    }



  delete_Text(&text);
  return 0;
}




void display_config(double * CONFIG, char * prob)
{
  printf("[%s] configurated:\n", prob);
  printf("  gamma     = %g\n", CONFIG[0]);

  printf("  MaxStp    = %d\n", (int)CONFIG[11]);
  printf("  TIME      = %g\n", CONFIG[12]);
  printf("  CFL       = %g\n", CONFIG[1]);
  printf("  eps       = %g\n", CONFIG[2]);
  printf("  tol       = %g\n", CONFIG[8]);
  printf("  scaling   = %g\n", CONFIG[19]);


  printf("  alpha1    = %g\tthe first  limiter of the slope\n", CONFIG[3]);
  printf("  alpha2    = %g\tthe second limiter of the slope\n", CONFIG[4]);
  printf("  beta1     = %g\tthe first  parameter of the monitor function\n", CONFIG[5]);
  printf("  beta2     = %g\tthe second parameter of the monitor function\n", CONFIG[6]);
  printf("  modfr     = %g\n", CONFIG[7]);
  printf("  threshold = %g\n", CONFIG[9]);
  printf("  thickness = %g\n", CONFIG[10]);

  if(CONFIG[13])
    CONFIG[13] = 1.0;
  else
    CONFIG[13] = 0.0;
  int adp = (int)CONFIG[13];
  if(adp)
    printf("  ADAPTIVE       mesh\n");
  else
    printf("  FIXED          mesh\n");

  int bod = (int)CONFIG[14];
  if(bod < 0)
    printf("  REFLECTION     boundary condition\n");
  else if(bod > 0)
    printf("  INFINITY       boundary condition\n");
  else
    printf("  PERIODIC       boundary condition\n");

  if(CONFIG[15])
    CONFIG[15] = 1.0;
  else
    CONFIG[15] = 0.0;
  int Primative = (int)CONFIG[15];
  if(Primative)
    printf("  PRIMITIVE       variables\n");
  else
    printf("  CONSERVATION    variables\n");

  if(CONFIG[16])
    CONFIG[16] = 1.0;
  else
    CONFIG[16] = 0.0;
  int Deri = (int)CONFIG[16];
  if(Deri)
    printf("  WENO-type     direvatives\n");
  else
    printf("  SMOOTH        direvatives\n");

  if(CONFIG[17])
    CONFIG[17] = 1.0;
  else
    CONFIG[17] = 0.0;
  int Limiter = (int)CONFIG[17];
  if(Limiter)
    printf("  Limiter              used\n");
  else
    printf("  Limiter          NOT used\n");

  if(CONFIG[18])
    CONFIG[18] = 1.0;
  else
    CONFIG[18] = 0.0;
  int Decomp = (int)CONFIG[18];
  if(Decomp)
    printf("  Decomposition        used\n");
  else
    printf("  Decomposition    NOT used\n");
  printf("\n");
}



int initialize(int nInitValue,  realArray InitValue[nInitValue], int sizeInitValue[nInitValue], char * err_msg, char addInitValue[nInitValue][L_STR], char * prob)
{
  int err_code = 02;
  FILE * fp_data;
  char extend[5] = ".txt\0";
  char add[L_STR+L_STR];
  strcpy(add, "../DATA/\0");
  int len_prob = 0, len_IV = 0, i, j, k, l;
  int it, idx;

  int n_digit_real = 15;
  char digit_real[] = {'0','1','2','3','4','5','6','7','8','9','-','+','.','e','E'};
  int n_space = 4;
  char space[] = {' ', '\t', '\n', '\r'};
  int read_state;
  char buffer[L_BUFF], number[L_BUFF], end_mark;
  int sign;

  while(prob[len_prob] != '\0')
    ++len_prob;
  for(it = 0; it < nInitValue; ++it)
  {
    len_IV = 0;
    while(addInitValue[it][len_IV] != '\0')
      ++len_IV;


    idx = 8;
    for(i = 0; i < len_prob; ++i)
      add[idx++] = prob[i];
    add[idx++] = '/';
    for(i = 0; i < len_IV; ++i)
      add[idx++] = addInitValue[it][i];
    for(i = 0; i < len_prob; ++i)
      add[idx++] = prob[i];
    for(i = 0; i < 5; ++i)
      add[idx++] = extend[i];




    if((fp_data = fopen(add, "r")) == 0)
    {
      sprintf(err_msg, "Cannot open initial value data file: %s!\n", add);
      return err_code;
    }

    init_realArray(InitValue +it);
    read_state = real_read(InitValue +it, &end_mark, fp_data, buffer, number, '#', EOF, space, n_space, digit_real, n_digit_real, err_msg);
    if(read_state)
    {
      //strcat(err_msg, sprintf(" (the %d-th item)\n", read_state));
      sprintf(err_msg, "%s (reading %s's %d-th item)\n", err_msg, addInitValue[it], read_state);
      for(k = 0; k <= it; ++k)
	delete_realArray(InitValue +k);
      return err_code + ((it+1) << bit_shift);
    }

    //display_realArray(InitValue +it);
    //printf("\n");

    //delete_realArray(InitValue +it);
    fclose(fp_data);
    sizeInitValue[it] = (InitValue[it].n_box-1)*InitValue[it].box_size+InitValue[it].tail_capacity;
  }

  return 0;
}
