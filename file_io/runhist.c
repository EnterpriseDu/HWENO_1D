#include <stdio.h>
#include <stdlib.h>

#include "file_io.h"


void init_runList(runList * runhist)
{
  runhist->length = 0;
  runhist->head = NULL;
  runhist->current = NULL;
  runhist->tail = NULL;
}

void insert_runList(runList * runhist)
{
  int j;

  if(runhist->length)
  {
    runhist->tail->next = (runNode *)malloc(sizeof(runNode));
    runhist->tail = runhist->tail->next;
  }
  else
  {
    runhist->tail = (runNode *)malloc(sizeof(runNode));
    runhist->head = runhist->tail;
  }
  ++(runhist->length);

  runhist->tail->trouble0 = NULL;
  runhist->tail->trouble1 = NULL;
  runhist->tail->next = NULL;
  for(j = 0; j < 6; ++j)
    runhist->tail->RcstrState[j] = '#';
  for(j = 0; j < 12; ++j)
    runhist->tail->RcstrErr[j] = 0.0;
  runhist->tail->time[0] = 0.0;
  runhist->tail->time[1] = 0.0;
}

void locate_runList(int p, runList * runhist)
{
  int j;
  runNode * point;

  point = runhist->head;
  for(j = 1; j < p; ++j)
    if(point->next)
      point = point->next;

  runhist->current = point;
}

void delete_runList(runList * runhist)
{
  int j, L = runhist->length;

  for(j = 0; j < L; ++j)
  {
    locate_runList(runhist->length-1, runhist);

    //printf("%c\n", runhist->tail->RcstrState[0]);
    if(runhist->tail->trouble0)
      free(runhist->tail->trouble0);
    if(runhist->tail->trouble1)
      free(runhist->tail->trouble1);
    free(runhist->tail);
    runhist->tail = runhist->current;
    --(runhist->length);

    if(runhist->length)
      runhist->tail->next = NULL;
    else
      runhist->head = NULL;
  }

  runhist->tail = runhist->head;
  locate_runList(0, runhist);
}


int check_runList(runList * runhist)
{
  int count = 0;
  int k = 0;
  
  runhist ->current = runhist->head;
  
  while(runhist->current)
  {
    if(runhist->tail == runhist->current)
      k = count+1;
    
    ++count;
    runhist->current = runhist->current->next;
  }

  if(runhist->length - k)
    return runhist->length - k;
  
  return runhist->length - count;
}
