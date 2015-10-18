#include <stdio.h>
#include <stdlib.h>


void matrix_distruction(int m, int n, double * pointer[])
{

  int j = 0;

  for(j = 0; j < m; ++j)
  {
    if(!pointer[j])
    {
      free(pointer[j]);
      pointer[j] = NULL;
    }
  }

}
