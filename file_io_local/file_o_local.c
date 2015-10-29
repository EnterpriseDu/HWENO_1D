#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "file_io.h"


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
