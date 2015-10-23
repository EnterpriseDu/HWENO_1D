#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "file_io.h"




extern double * RHO0;
extern double * U0;
extern double * V0;
extern double * P0;
extern double * X0;
extern double * Y0;



/*
 * vvM = 0 or vvM = K
 *
 *
 *
 */
void file_write_data
(int m, int start, int vvM, double * data[], char * source, char * sol_name)
{
  FILE * fp_write;
  char add_data[100] = "";

//===================Write solution File=========================

  strcat(add_data, "../SOLUTION/");
  strcat(add_data, sol_name);
  strcat(add_data, "/");
  strcat(add_data, source);
  strcat(add_data, "_0000.txt");

  int p4, p3, p2, p1, ptmp, len;

  int j = 0, i = 0, k = 0;
  for(k = start; k <= vvM; ++k)
  {
    p4 = k%10;
    ptmp = (int)(k / 10);
    p3 = ptmp % 10;
    ptmp = (int)(ptmp / 10);
    p2 = ptmp % 10;
    p1 = (int)(ptmp / 10);

    len = strlen(add_data);
    add_data[len-5] = (char)(p4+48);
    add_data[len-6] = (char)(p3+48);
    add_data[len-7] = (char)(p2+48);
    add_data[len-8] = (char)(p1+48);

    if((fp_write = fopen(add_data, "w")) == 0)
    {
      printf("Cannot open solution output file: %s!\n", add_data);
      exit(1);
    }
    for(j = 0; j < m; ++j)
    {
      fprintf(fp_write, "%.18lf\t", data[k][j]);
    }
    fclose(fp_write);
  }
}




/*
 *
 *
 *
 *
 */
void file_write_log
(int m, int n, int K, double scaling, double * CONFIG, double * OPT, runList * runhist, char * scheme, char * prob, char * sol_name)
{
  FILE * fp_write;
  char add_log[100] = "";

  time_t t;
  struct tm * local_time;

  double sum_tau = 0.0, sum_cpu = 0.0;

  int j = 0, i = 0, k = 0;



  if(check_runList(runhist))
  {
    printf("In the function [file_write_log]: the runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runList(runhist));
    exit(100);
  }
  if(K - runhist->length)
  {
    printf("In the function [file_write_log]: after %d steps of computation, the number of the runNodes is %d.\n\n", K, runhist->length);
    exit(100);
  }

  

  strcat(add_log, "../SOLUTION/");
  strcat(add_log, sol_name);
  strcat(add_log, "/log.txt");

  if((fp_write = fopen(add_log, "w")) == 0)
  {
    printf("Cannot open log output file: %s!\n", add_log);
    exit(1);
  }

  fprintf(fp_write, "[%s | %s] initialized with %d x-grids, %d y-grids.\n", scheme, prob, m, n);

  t=time(NULL);
  local_time=localtime(&t);
  fprintf(fp_write, "The computation finished at [%4d.%02d.%02d|%02d:%02d:%02d].\n", local_time->tm_year+1900, local_time->tm_mon+1, local_time->tm_mday, local_time->tm_hour, local_time->tm_min, local_time->tm_sec);

  
  runhist->current = runhist->head;
  while(runhist->current)
  {
    sum_cpu += runhist->current->time[1];
    runhist->current = runhist->current->next;
  }
  fprintf(fp_write, "%d steps computed and %g seconds CPU time used.\n\n", K, sum_cpu);
  sum_cpu = 0.0;

  fprintf(fp_write, "  gamma  = %g\n", CONFIG[0]);
  fprintf(fp_write, "  CFL    = %g\n", CONFIG[1]);
  fprintf(fp_write, "  eps    = %g\n", CONFIG[2]);
  fprintf(fp_write, "  alpha1 = %g\tthe first  limiter of the slope\n", CONFIG[3]);
  fprintf(fp_write, "  alpha2 = %g\tthe second limiter of the slope\n", CONFIG[4]);
  fprintf(fp_write, "  beta1  = %g\tthe first  parameter of the monitor function\n", CONFIG[5]);
  fprintf(fp_write, "  beta2  = %g\tthe second parameter of the monitor function\n", CONFIG[6]);
  fprintf(fp_write, "  modfr  = %g\n", CONFIG[7]);
  fprintf(fp_write, "  tol    = %g\n", CONFIG[8]);

  fprintf(fp_write, "  MaxStp = %d\n", (int)OPT[0]);
  fprintf(fp_write, "  TIME   = %g\n", OPT[1]);
  int inter_data = (int)OPT[2];
  if(inter_data)
    fprintf(fp_write, "  inter-data     kept\n");
  else
    fprintf(fp_write, "  inter-data NOT kept\n");

  int adp = (int)OPT[3];
  if(adp)
    fprintf(fp_write, "  ADAPTIVE       mesh\n");
  else
    fprintf(fp_write, "  FIXED          mesh\n");

  int bod = (int)OPT[4];
  if(bod < 0)
    fprintf(fp_write, "  REFLECTION     boundary condition\n");
  else if(bod > 0)
    fprintf(fp_write, "  INFINITY       boundary condition\n");
  else
    fprintf(fp_write, "  PERIODIC       boundary condition\n");

  int Riemann = (int)OPT[5];
  if(Riemann)
    fprintf(fp_write, "  PRIMITIVE      variables\n");
  else
    fprintf(fp_write, "  CONSERVATION   variables\n");

  int WENOD = (int)OPT[6];
  if(WENOD)
    fprintf(fp_write, "  WENO-type      direvatives\n");
  else
    fprintf(fp_write, "  SMOOTH         direvatives\n");

  int limiter = (int)OPT[7];
  if(limiter)
    fprintf(fp_write, "  LIMITER        used\n");
  else
    fprintf(fp_write, "  limiter    NOT used\n");
  fprintf(fp_write, "\n");

  
  k = 0;
  runhist->current = runhist->head;
  while(runhist->current)
  {
    ++k;
    fprintf(fp_write, "\n-------------%04d-------------\n", k);
    runhist->current->time[0] = runhist->current->time[0]/scaling;
    sum_tau += runhist->current->time[0];
    sum_cpu += runhist->current->time[1];
    fprintf(fp_write, " time step length\n");
    fprintf(fp_write, "    t^n=%3.12lf      tau=%3.12lf\n", sum_tau, runhist->current->time[0]);
    fprintf(fp_write, "  cpu time:\n");
    fprintf(fp_write, "    total=%lf      local=%lf\n", sum_cpu, runhist->current->time[1]);

    if(adp)
    {
      fprintf(fp_write, "  redistribution:\n");
      fprintf(fp_write, "    %c  %6d  %g | ", runhist->current->RcstrState[0], (int)runhist->current->RcstrErr[0], runhist->current->RcstrErr[1]);
      fprintf(fp_write, "%c  %6d  %g | ",     runhist->current->RcstrState[1], (int)runhist->current->RcstrErr[2], runhist->current->RcstrErr[3]);
      fprintf(fp_write, "%c  %6d  %g | ",     runhist->current->RcstrState[2], (int)runhist->current->RcstrErr[4], runhist->current->RcstrErr[5]);
      fprintf(fp_write, "%c  %6d  %g\n",      runhist->current->RcstrState[3], (int)runhist->current->RcstrErr[6], runhist->current->RcstrErr[7]);
      fprintf(fp_write, "    %c  %6d  %g | ", runhist->current->RcstrState[4], (int)runhist->current->RcstrErr[8], runhist->current->RcstrErr[9]);
      fprintf(fp_write, "%c  %6d  %g\n",      runhist->current->RcstrState[5], (int)runhist->current->RcstrErr[10], runhist->current->RcstrErr[11]);
    }
    runhist->current = runhist->current->next;
    if(k > K)
      break;
  }
  
  fclose(fp_write);
}


void file_write_trouble(int m, int K, runList * runhist, char * sol_name)
{
  FILE * fp_write;
  char add_data[100] = "";
  int j = 0, i = 0, k = 0;

  


  if(check_runList(runhist))
  {
    printf("In the function [file_write_trouble]: the runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runList(runhist));
    exit(100);
  }
  if(K - runhist->length)
  {
    printf("In the function [file_write_trouble]: after %d steps of computation, the number of the runNodes is %d.\n\n", K, runhist->length);
    exit(100);
  }

//===================Write solution File=========================

  strcat(add_data, "../SOLUTION/");
  strcat(add_data, sol_name);
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







void write_column(int m, double data[], char * source, char * sol_name)
{
  FILE * fp_write;
  char add_data[100] = "";

//===================Write solution File=========================

  strcat(add_data, "../SOLUTION/");
  strcat(add_data, sol_name);
  strcat(add_data, "/");
  strcat(add_data, source);
  strcat(add_data, ".txt");


  int j = 0, i = 0;

  if((fp_write = fopen(add_data, "w")) == 0)
  {
    printf("Cannot open solution output file: %s!\n", add_data);
    exit(1);
  }

  for(j = 0; j < m; ++j)
    fprintf(fp_write, "%.10lf\n", data[j]);


  fclose(fp_write);
}
