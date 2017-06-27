#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "file_io.h"
#ifndef L_STR
#include "file_io_local.h"
#endif



int DATA_OUTPUT(char * err_msg, int nDATA, char addDATA[nDATA][L_STR], int sizeDATA[nDATA], double * DATA[nDATA], int idx[nDATA], int flagDATA[nDATA], char * directory)
{
  int err_code = 20;
  FILE * fp_write;
  char add_data[L_STR+L_STR] = "";
  char strIDX[L_STR];
  int j = 0, i = 0, k = 0, it;

//===================Write solution File=========================

  for(it = 0; it < nDATA; ++it)
  {
    if(!flagDATA[it])
      continue;
    strcpy(add_data, directory);
    //strcat(add_data, "/");
    strcat(add_data, addDATA[it]);
    sprintf(strIDX, "_%04d", idx[it]);
    strcat(add_data, strIDX);
    strcat(add_data, ".txt\0");

    if((fp_write = fopen(add_data, "w")) == 0)
    {
      sprintf(err_msg, "Cannot open solution output file: %s!\n", add_data);
      return err_code + (it << bit_shift);
    }
    for(j = 0; j < sizeDATA[it]; ++j)
      fprintf(fp_write, "%.18lf\t", DATA[it][j]);
    fclose(fp_write);
  }

  return 0;
}


int LOG_OUTPUT(char * err_msg, runHist * runhist, char ITEM[N_CONF][L_STR], int already_read[N_CONF], double CONFIG[N_CONF], int m, int n, int K, char * scheme, char * version, char * prob, char * directory)
{
  int err_code = 30;
  FILE * fp_write;
  char add_log[L_STR+L_STR] = "";
  int j = 0, i = 0, k = 0, it, len, LEN;
  time_t t;
  struct tm * local_time;
  double sum_tau = 0.0, sum_cpu = 0.0;
  int flag_time[4] = {0,0,0,0}, flag_extra[2] = {0,0};

  if(check_runHist(runhist))
  {
    printf("In the function [file_write_log]: the runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runHist(runhist));
    return err_code + (10 << bit_shift);
  }
  if(K - runhist->length)
  {
    printf("In the function [file_write_log]: after %d steps of computation, the number of the runNodes is %d.\n\n", K, runhist->length);
    return err_code + (20 << bit_shift);
  }
//===================Write log File=========================


    strcpy(add_log, directory);
    strcat(add_log, "log.txt\0");


    if((fp_write = fopen(add_log, "w")) == 0)
    {
      sprintf(err_msg, "Cannot open solution output file: %s!\n", add_log);
      return err_code + (0 << bit_shift);
    }

    fprintf(fp_write, "[%s | %s] initialized with %d x-grids, %d y-grids.\n", scheme, prob, m, n);
    fprintf(fp_write, "The present version is [%s].\n", version);

    t=time(NULL);
    local_time=localtime(&t);
    fprintf(fp_write, "The computation finished at [%4d.%02d.%02d|%02d:%02d:%02d].\n", local_time->tm_year+1900, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);

    LEN = 0;
    for(it = 0; it < N_CONF; ++it)
      if(LEN < strlen(ITEM[it]))
	LEN = strlen(ITEM[it]);
    for(it = 0; it < N_CONF; ++it)
    {
      len = fprintf(fp_write, "%s", ITEM[it]);
      for(j = 0; j < LEN+5-len; ++j)
	fprintf(fp_write, " ");
      fprintf(fp_write, "%g", CONFIG[it]);
      if(already_read[it] < 0)
	fprintf(fp_write, "  (default)");
      fprintf(fp_write, "\n");
    }

    //flag_time = {0,0,0,0}; //tau, T, current cpu, total cpu
    //flag_extra = {0,0}; //int, double
    sum_cpu = write_runHist(runhist, fp_write, 0, flag_time, flag_extra, (int)(CONFIG[13]), CONFIG[19], '\t');
    fprintf(fp_write, "The cost of cpu time is %.18lf.\n\n", sum_cpu);
    //flag_time[1] = 1; // output T
    //flag_extra[0] = 1; // output the INT
    //sum_cpu = write_runHist(runhist, fp_write, 0, flag_time, flag_extra, (int)(CONFIG[13]), CONFIG[19], '\t');

    fclose(fp_write);

  return 0;
}


/*
void file_write_trouble(int m, int K, runHist *runhist, char * directory)
{
  FILE * fp_write;
  char add_data[100] = "";
  int j = 0, i = 0, k = 0;

  
  if(check_runHist(runhist))
  {
    printf("In the function [file_write_trouble]: the runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runHist(runhist));
    exit(100);
  }
  if(K - runhist->length)
  {
    printf("In the function [file_write_trouble]: after %d steps of computation, the number of the runNodes is %d.\n\n", K, runhist->length);
    exit(100);
  }

//===================Write solution File=========================

  strcat(add_data, directory);
  strcat(add_data, "/trouble0");
  strcat(add_data, ".txt");

  if((fp_write = fopen(add_data, "w")) == 0)
  {
    printf("Cannot open solution output file: %s!\n", add_data);
    exit(1);
  }
  
  k = 0;
  runhist->current = runhist->head;
  while(runhist->current)
  {
    ++k;

    if(runhist->current->trouble0)
    {
      for(j = 0; j < m; ++j)
	fprintf(fp_write, "%d\t", runhist->current->trouble0[j]);
      fprintf(fp_write, "\n");
    }
    
    runhist->current = runhist->current->next;
    if(k > K)
      break;
  }



  k = 0;
  while(add_data[k++]);
  //++k;
  add_data[k-6] = '1';

  if((fp_write = fopen(add_data, "w")) == 0)
  {
    printf("Cannot open solution output file: %s!\n", add_data);
    exit(1);
  }
  
  k = 0;
  runhist->current = runhist->head;
  while(runhist->current)
  {
    ++k;

    if(runhist->current->trouble1)
    {
      for(j = 0; j < m; ++j)
	fprintf(fp_write, "%d\t", runhist->current->trouble1[j]);
      fprintf(fp_write, "\n");
    }
    
    runhist->current = runhist->current->next;
    if(k > K)
      break;
  }
  
}
//*/




