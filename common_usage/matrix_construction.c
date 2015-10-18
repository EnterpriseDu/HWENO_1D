#include <stdio.h>
#include <stdlib.h>


int matrix_construction(int m, int n, double * pointer[])
{

  int j = 0;

  for(j = 0; j < m; ++j)
  {
    if(!pointer[j])
    {
      free(pointer[j]);
      pointer[j] = NULL;
    }

    pointer[j] = (double *)malloc(sizeof(double) * n);

    if(!pointer[j])
      return 0;
  }

  return 1;
}
