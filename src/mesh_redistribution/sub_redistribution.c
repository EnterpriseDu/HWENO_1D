#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mesh_redistribution.h"


/*
 * This function redistribute the 2-D quadrangular mesh
 * using the method of Winslow's invariation.
 */
void sub_redistribution
(int m, int n, double * x[m+1], double * y[m+1], double * OMEGA[m],
 int max_it, double tol, double eps, double * number, char * charact)
{
  int const DBG = 0; //compile the debug code or not

  int i = 0, j = 0, k = 0, idx, M = max_it;
  double err;

  //variables used in the CG algorithm
  double * res, * res_tilde, * p, * Ap, step_a, step_b, inner_prod;

  res_tilde = (double *)malloc(sizeof(double) * (m-1));
  res = (double *)malloc(sizeof(double) * (m-1));
  p   = (double *)malloc(sizeof(double) * (m-1));
  Ap  = (double *)malloc(sizeof(double) * (m-1));

  for(j = 1; j < m; ++j)
  {
    res[j-1] = (OMEGA[j-1][0]+OMEGA[j][0])*x[j][0];
    res[j-1] -= OMEGA[j-1][0]*x[j-1][0];
    res[j-1] -= OMEGA[j][0]*x[j+1][0];
    p[j-1] = res[j-1];
  }
  err = 0.0;
  for(j = 0; j < m-1; ++j)
    err += res[j]*res[j];
  //err = sqrt(err) / (double)m;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute the bottom boundary with initial error=%lf.\n", err);
    number[2] = 0.0;
    number[3] = 0.0;
    charact[0] = 'O';
  }

  for(k = 1; k < M; ++k)
  {
    for(j = 1; j < m; ++j)
    {
      Ap[j-1] = -(OMEGA[j-1][0]+OMEGA[j][0])*p[j-1];
      if(j != 1)
	Ap[j-1] += OMEGA[j-1][0]*p[j-2];
      if(j != m-1)
	Ap[j-1] += OMEGA[j][0]*p[j];
    }

    step_a = 0.0;
    inner_prod = 0.0;
    for(j = 0; j < m-1; ++j)
    {
      step_a += res[j]*res[j];
      inner_prod += p[j]*Ap[j];
    }
    step_a = step_a/inner_prod;

    for(j = 1; j < m; ++j)
	x[j][0] += step_a * p[j-1];

    for(j = 0; j < m-1; ++j)
      res_tilde[j] = res[j] - step_a * Ap[j];

    step_b = 0.0;
    inner_prod = 0.0;
    for(j = 0; j < m-1; ++j)
    {
      step_b += res_tilde[j]*res_tilde[j];
      inner_prod += res[j]*res[j];
    }
    step_b = step_b/inner_prod;

    for(j = 0; j < m-1; ++j)
    {
      res[j] = res_tilde[j];
      p[j] = res[j] + step_b * p[j];
    }

    err = 0.0;
    for(j = 0; j < m-1; ++j)
       err += res[j]*res[j];

    //err = sqrt(err) / (double)m;
    //printf("err: %g\n", err);
    if(err < tol)
      break;
  }

  if(M)
  {
    number[2] = (double)k;
    number[3] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for the bottom boundary complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[0] = '$';
    }
    else
    {
      printf("Maximum number of iteration for the bottom boundary reached!\nThe total error is still %g\n", err);
      charact[0] = '#';
    }
  }


  M = max_it;
  for(j = 1; j < m; ++j)
  {
    res[j-1] = (OMEGA[j-1][n-1]+OMEGA[j][n-1])*x[j][n];
    res[j-1] -= OMEGA[j-1][n-1]*x[j-1][n];
    res[j-1] -= OMEGA[j][n-1]*x[j+1][n];
    p[j-1] = res[j-1];
  }
  err = 0.0;
  for(j = 0; j < m-1; ++j)
    err += res[j]*res[j];
  //err = sqrt(err) / (double)m;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute the top boundary with initial error=%lf.\n", err);
    number[4] = 0.0;
    number[5] = 0.0;
    charact[1] = 'O';
  }

  for(k = 1; k < M; ++k)
  {
    for(j = 1; j < m; ++j)
    {
      Ap[j-1] = -(OMEGA[j-1][n-1]+OMEGA[j][n-1])*p[j-1];
      if(j != 1)
	Ap[j-1] += OMEGA[j-1][n-1]*p[j-2];
      if(j != m-1)
	Ap[j-1] += OMEGA[j][n-1]*p[j];
    }

    step_a = 0.0;
    inner_prod = 0.0;
    for(j = 0; j < m-1; ++j)
    {
      step_a += res[j]*res[j];
      inner_prod += p[j]*Ap[j];
    }
    step_a = step_a/inner_prod;

    for(j = 1; j < m; ++j)
	x[j][n] += step_a * p[j-1];

    for(j = 0; j < m-1; ++j)
      res_tilde[j] = res[j] - step_a * Ap[j];

    step_b = 0.0;
    inner_prod = 0.0;
    for(j = 0; j < m-1; ++j)
    {
      step_b += res_tilde[j]*res_tilde[j];
      inner_prod += res[j]*res[j];
    }
    step_b = step_b/inner_prod;

    for(j = 0; j < m-1; ++j)
    {
      res[j] = res_tilde[j];
      p[j] = res[j] + step_b * p[j];
    }

    err = 0.0;
    for(j = 0; j < m-1; ++j)
       err += res[j]*res[j];

    //err = sqrt(err) / (double)m;
    //printf("err: %g\n", err);
    if(err < tol)
      break;
  }
  
  if(M)
  {
    number[4] = (double)k;
    number[5] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for the top boundary complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[1] = '$';
    }
    else
    {
      printf("Maximum number of iteration for the top boundary reached!\nThe total error is still %g\n", err);
      charact[1] = '#';
    }
  }




  free(res);
  free(res_tilde);
  free(p);
  free(Ap);
  res_tilde = (double *)malloc(sizeof(double) * (n-1));
  res = (double *)malloc(sizeof(double) * (n-1));
  p   = (double *)malloc(sizeof(double) * (n-1));
  Ap  = (double *)malloc(sizeof(double) * (n-1));


  M = max_it;
  for(i = 1; i < n; ++i)
  {
    res[i-1] = (OMEGA[0][i-1]+OMEGA[0][i])*y[0][i];
    res[i-1] -= OMEGA[0][i-1]*y[0][i-1];
    res[i-1] -= OMEGA[0][i]*y[0][i+1];
    p[i-1] = res[i-1];
  }
  err = 0.0;
  for(j = 0; j < n-1; ++j)
    err += res[j]*res[j];
  //err = sqrt(err) / (double)n;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute the left boundary with initial error=%lf.\n", err);
    number[6] = 0.0;
    number[7] = 0.0;
    charact[2] = 'O';
  }

  for(k = 1; k < M; ++k)
  {
    for(i = 1; i < n; ++i)
    {
      Ap[i-1] = -(OMEGA[0][i-1]+OMEGA[0][i])*p[i-1];
      if(i != 1)
	Ap[i-1] += OMEGA[0][i-1]*p[i-2];
      if(i != n-1)
	Ap[i-1] += OMEGA[0][i]*p[i];
    }

    step_a = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < n-1; ++i)
    {
      step_a += res[i]*res[i];
      inner_prod += p[i]*Ap[i];
    }
    step_a = step_a/inner_prod;

    for(i = 1; i < n; ++i)
	y[0][i] += step_a * p[i-1];

    for(i = 0; i < n-1; ++i)
      res_tilde[i] = res[i] - step_a * Ap[i];

    step_b = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < n-1; ++i)
    {
      step_b += res_tilde[i]*res_tilde[i];
      inner_prod += res[i]*res[i];
    }
    step_b = step_b/inner_prod;

    for(i = 0; i < n-1; ++i)
    {
      res[i] = res_tilde[i];
      p[i] = res[i] + step_b * p[i];
    }

    err = 0.0;
    for(i = 0; i < n-1; ++i)
       err += res[i]*res[i];

    //err = sqrt(err) / (double)n;
    //printf("err: %g\n", err);
    if(err < tol)
      break;
  }

  if(M)
  {
    number[6] = (double)k;
    number[7] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for the left boundary complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[2] = '$';
    }
    else
    {
      printf("Maximum number of iteration for the left boundary reached!\nThe total error is still %g\n", err);
      charact[2] = '#';
    }
  }



  M = max_it;
  for(i = 1; i < n; ++i)
  {
    res[i-1] = (OMEGA[m-1][i-1]+OMEGA[m-1][i])*y[m][i];
    res[i-1] -= OMEGA[m-1][i-1]*y[m][i-1];
    res[i-1] -= OMEGA[m-1][i]*y[m][i+1];
    p[i-1] = res[i-1];
  }
  err = 0.0;
  for(j = 0; j < n-1; ++j)
    err += res[j]*res[j];
  //err = sqrt(err) / (double)n;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute the right boundary with initial error=%lf.\n", err);
    number[8] = 0.0;
    number[9] = 0.0;
    charact[3] = 'O';
  }

  for(k = 1; k < M; ++k)
  {
    for(i = 1; i < n; ++i)
    {
      Ap[i-1] = -(OMEGA[m-1][i-1]+OMEGA[m-1][i])*p[i-1];
      if(i != 1)
	Ap[i-1] += OMEGA[m-1][i-1]*p[i-2];
      if(i != n-1)
	Ap[i-1] += OMEGA[m-1][i]*p[i];
    }

    step_a = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < n-1; ++i)
    {
      step_a += res[i]*res[i];
      inner_prod += p[i]*Ap[i];
    }
    step_a = step_a/inner_prod;

    for(i = 1; i < n; ++i)
	y[m][i] += step_a * p[i-1];

    for(i = 0; i < n-1; ++i)
      res_tilde[i] = res[i] - step_a * Ap[i];

    step_b = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < n-1; ++i)
    {
      step_b += res_tilde[i]*res_tilde[i];
      inner_prod += res[i]*res[i];
    }
    step_b = step_b/inner_prod;

    for(i = 0; i < n-1; ++i)
    {
      res[i] = res_tilde[i];
      p[i] = res[i] + step_b * p[i];
    }

    err = 0.0;
    for(i = 0; i < n-1; ++i)
       err += res[i]*res[i];

    //err = sqrt(err) / (double)n;
    //printf("err: %g\n", err);
    if(err < tol)
      break;
  }

  if(M)
  {
    number[8] = (double)k;
    number[9] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for the right boundary complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[3] = '$';
    }
    else
    {
      printf("Maximum number of iteration for the right boundary reached!\nThe total error is still %g\n", err);
      charact[3] = '#';
    }
  }

  free(res);
  free(res_tilde);
  free(p);
  free(Ap);
}
