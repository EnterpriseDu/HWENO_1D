#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mesh_redistribution.h"


/*
 * This function redistribute the 2-D quadrangular mesh
 * by using Winslow's variation method.
 *
 * As is well known, Winslow's variation method request
 * to solve a system of decoupled eliptic equations:
 *
 *     (omega * x_xi)_xi + (omega * x_eta)_eta = 0
 *     (omega * y_xi)_xi + (omega * y_eta)_eta = 0,
 *
 * where {xi, eta} is the logic coordinates. Each of the
 *  above equations can be discretized into a linear system
 * of equations with a five-diagnal coefficient matrix A
 * by five-point difference method.
 *
 * In this function, we use pre-conditioner conjugate
 * gradient (PCG) method, with the pre-conditioner B to
 * be a diagnal matrix whose diagnal entries is the inverse
 * of the diagnal entries of A.
 */
int redistribution
(int m, int n, double * x[m+1], double * y[m+1], double * OMEGA[m],
 // double * s_rho[m], double * s_u[m], double * s_v[m], double s_p[m],
 // double * t_rho[m], double * t_u[m], double * t_v[m], double t_p[m],
 int max_it, double tol, double eps, double * number, char * charact)
{
  int const DBG = 0; //compile the debug code or not
  int M = max_it;
  int i = 0, j = 0, k = 0, idx, len = (m-1)*(n-1);
  double err;

  //variables use in the CG algorithm
  double * res, * res_tilde, * p, * Ap, step_a, step_b, inner_prod;
  double * xv;  //rearrange x or y into a vector

  res_tilde = (double *)malloc(sizeof(double) * len);
  res = (double *)malloc(sizeof(double) * len);
  p   = (double *)malloc(sizeof(double) * len);
  Ap  = (double *)malloc(sizeof(double) * len);
  xv = (double *)malloc(len*sizeof(double));


  for(j = 1; j < m; ++j)
    for(i = 1; i < n; ++i)
      xv[j*n -j -n +i] = -x[j][i];
  Axfun(m, n, res, OMEGA, xv);
  for(j = 1; j < m; ++j)
  {
    idx = (j-1) * (n-1);
    res[idx] += 0.5*(OMEGA[j][0]+OMEGA[j-1][0])/(OMEGA[j-1][0]+OMEGA[j-1][1]+OMEGA[j][0]+OMEGA[j][1]) * x[j][0];
    res[idx+n-2] += 0.5*(OMEGA[j][n-1]+OMEGA[j-1][n-1])/(OMEGA[j-1][n-2]+OMEGA[j-1][n-1]+OMEGA[j][n-2]+OMEGA[j][n-1]) * x[j][n];
  }
  idx = (m-2) * (n-1);
  for(i = 1; i < n; ++i)
  {
    res[i-1] += 0.5*(OMEGA[0][i]+OMEGA[0][i-1])/(OMEGA[0][i-1]+OMEGA[0][i]+OMEGA[1][i-1]+OMEGA[1][i]) * x[0][i];
    res[idx+i-1] += 0.5*(OMEGA[m-1][i]+OMEGA[m-1][i-1])/(OMEGA[m-2][i-1]+OMEGA[m-2][i]+OMEGA[m-1][i-1]+OMEGA[m-1][i]) * x[m][i];
  }

  err = 0.0;
  for(i = 0; i < len; ++i)
    err += res[i]*res[i];
  //err = sqrt(err) / len;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute y with initial error=%lf.\n", err);
    number[12] = 0.0;
    number[13] = 0.0;
    charact[5] = 'O';
  }

  for(i = 0; i < len; ++i)
    p[i] = res[i];

  for(k = 1; k < M; ++k)
  {
    Axfun(m, n, Ap, OMEGA, p);
    step_a = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < len; ++i)
    {
      step_a += res[i]*res[i];
      inner_prod += p[i]*Ap[i];
    }
    step_a = step_a/inner_prod;

    for(j = 1; j < m; ++j)
      for(i = 1; i < n; ++i)
	x[j][i] += step_a * p[j*n - j - n + i];

    for(i = 0; i < len; ++i)
      res_tilde[i] = res[i] - step_a * Ap[i];

    step_b = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < len; ++i)
    {
      step_b += res_tilde[i]*res_tilde[i];
      inner_prod += res[i]*res[i];
    }
    step_b = step_b/inner_prod;

    for(i = 0; i < len; ++i)
    {
      res[i] = res_tilde[i];
      p[i] = res[i] + step_b * p[i];
    }

    err = 0.0;
    for(i = 0; i < len; ++i)
       err += res[i]*res[i];

    //err = sqrt(err) / len;
    if(err < tol)
      break;
  }

  if(M)
  {
    number[12] = (double)k;
    number[13] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for y complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[5] = '$';
    }
    else
    {
      printf("Maximum number of iteration for y reached!\nThe total error is still %g\n", err);
      charact[5] = '#';
    }
  }


  /*
   *
   *              |                           |                           |
   *              |                           |                           |
   *              |                           |                           |
   *         [j-1][i+1]                    [j][i+1]                   [j+1][i+1]
   * --------x(j-1,i+1)--------------------x(j,i+1)--------------------x(j+1,i)----------
   *        [j*n-j-n-n+i]               [j*n-j-n+i+1]                 [j*n-j+i]
   *              |                           |                           |
   *              |                           |                           |
   *              |      OMEGA(j-0.5,i+0.5)   |   OMEGA(j+0.5,i+0.5)      |
   *              |            [j-1][i]       |           [j][i]          |
   *              |                           |                           |
   *         [j-1][i]                       [j][i]                     [j+1][i]
   * --------x(j-1,i)-----------------------x(j,i)---------------------x(j+1,i)----------
   *      [j*n-j-n-n+i-1]                [j*n-j-n+i]                 [j*n-j+i-1]
   *              |                           |                           |
   *              |                           |                           |
   *              |      OMEGA(j-0.5,i-0.5)   |   OMEGA(j+0.5,i-0.5)      |
   *              |            [j-1][i-1]     |           [j][i-1]        |
   *              |                           |                           |
   *        [j-1][i-1]                     [j][i-1]                   [j+1][i-1]
   * -------x(j-1,i-1)--------------------x(j,i-1)--------------------x(j+1,i-1)----------
   *      [j*n-j-n-n+i-2]               [j*n-j-n+i-1]                [j*n-j+i-2]
   *              |                           |                           |
   *              |                           |                           |
   *              |                           |                           |
   *
   *
   *       |                           |
   *       |                           |
   *       |                           |
   *   [0][i+1]                    [1][i+1]
   *   x(0,i+1)--------------------x(1,i+1)---------
   *       |                          [i]
   *       |                           |
   *       |                           |
   *       |      OMEGA(-0.5,i+0.5)    |
   *       |             [0][i]        |
   *       |                           |
   *    [0][i]                      [1][i]
   *    x(0,i)----------------------x(1,i)----------
   *       |                         [i-1]
   *       |                           |
   *       |                           |
   *       |      OMEGA(-0.5,i-0.5)    |
   *       |             [0][i-1]      |
   *       |                           |
   *   [0][i-1]                     [1][i-1]
   *   x(0,i-1)--------------------x(1,i-1)--------
   *       |                         [i-2]
   *       |                           |
   *       |                           |
   *       |                           |
   *
   *
   *              |                           |
   *              |                           |
   *              |                           |
   *          [m-1][i+1]                   [m][i+1]
   * ---------x(m-1,i+1)-------------------x(m,i+1)
   *        [(m-2)(n-1)+i]                    |
   *              |                           |
   *              |                           |
   *              |      OMEGA(m-0.5,i+0.5)   |
   *              |            [m-1][i]       |
   *              |                           |
   *          [m-1][i]                      [m][i]
   * ---------x(m-1,i)----------------------x(m,i)
   *     [(m-2)(n-1)+i-1]                     |
   *              |                           |
   *              |                           |
   *              |      OMEGA(m-0.5,i-0.5)   |
   *              |            [m-1][i-1]     |
   *              |                           |
   *         [m-1][i-1]                    [m][i-1]
   * --------x(m-1,i-1)--------------------x(m,i-1)
   *      [(m-2)(n-1)+i-2]                    |
   *              |                           |
   *              |                           |
   *              |                           |
   *
   *
   *
   *              |                           |                           |
   *              |                           |                           |
   *              |                           |                           |
   *         [j-1][i]                       [j][i]                     [j+1][i]
   * --------x(j-1,i)-----------------------x(j,i)---------------------x(j+1,i)----------
   *      [j*n-j-n-n+i-1]                [j*n-j-n+i]                 [j*n-j+i-1]
   *              |                           |                           |
   *              |                           |                           |
   *              |      OMEGA(j-0.5,0.5)   |   OMEGA(j+0.5,0.5)      |
   *              |            [j-1][0]     |           [j][0]        |
   *              |                           |                           |
   *        [j-1][0]                     [j][0]                   [j+1][i-1]
   * -------x(j-1,0)--------------------x(j,0)--------------------x(j+1,i-1)----------
   *      [j*n-j-n-n+i-2]               [j*n-j-n+i-1]                [j*n-j+i-2]
   *
   *
   *              |                           |                           |
   *              |                           |                           |
   *              |                           |                           |
   *         [j-1][i+1]                    [j][i+1]                   [j+1][i+1]
   * --------x(j-1,i+1)--------------------x(j,i+1)--------------------x(j+1,i)----------
   *        [j*n-j-n-n+i]               [j*n-j-n+i+1]                 [j*n-j+i]
   *              |                           |                           |
   *              |                           |                           |
   *              |      OMEGA(j-0.5,i+0.5)   |   OMEGA(j+0.5,i+0.5)      |
   *              |            [j-1][i]       |           [j][i]          |
   *              |                           |                           |
   *         [j-1][i]                       [j][i]                     [j+1][i]
   * --------x(j-1,i)-----------------------x(j,i)---------------------x(j+1,i)----------
   *      [j*n-j-n-n+i-1]                [j*n-j-n+i]                 [j*n-j+i-1]
   *              |                           |                           |
   *              |                           |                           |
   *              |      OMEGA(j-0.5,i-0.5)   |   OMEGA(j+0.5,i-0.5)      |
   *              |            [j-1][i-1]     |           [j][i-1]        |
   *              |                           |                           |
   *        [j-1][i-1]                     [j][i-1]                   [j+1][i-1]
   * -------x(j-1,i-1)--------------------x(j,i-1)--------------------x(j+1,i-1)----------
   *      [j*n-j-n-n+i-2]               [j*n-j-n+i-1]                [j*n-j+i-2]
   *              |                           |                           |
   *              |                           |                           |
   *              |                           |                           |
   */





  for(i = 1; i < n; ++i)
    for(j = 1; j < m; ++j)
      xv[i*m -i -m +j] = -y[j][i];
  Ayfun(m, n, res, OMEGA, xv);
  for(i = 1; i < n; ++i)
  {
    idx = (i-1) * (m-1);
    res[idx] += 0.5*(OMEGA[0][i]+OMEGA[0][i-1])/(OMEGA[0][i-1]+OMEGA[1][i-1]+OMEGA[0][i]+OMEGA[1][i]) * y[0][i];
    res[idx+m-2] += 0.5*(OMEGA[m-1][i]+OMEGA[m-1][i-1])/(OMEGA[m-2][i-1]+OMEGA[m-1][i-1]+OMEGA[m-2][i]+OMEGA[m-1][i]) * y[m][i];
  }
  idx = (n-2) * (m-1);
  for(j = 1; j < m; ++j)
  {
    res[j-1] += 0.5*(OMEGA[j][0]+OMEGA[j-1][0])/(OMEGA[j-1][0]+OMEGA[j][0]+OMEGA[j-1][1]+OMEGA[j][1]) * y[j][0];
    res[idx+j-1] += 0.5*(OMEGA[j][n-1]+OMEGA[j-1][n-1])/(OMEGA[j-1][n-2]+OMEGA[j][n-2]+OMEGA[j-1][n-1]+OMEGA[j][n-1]) * y[j][n];
  }

  M = max_it;
  err = 0.0;
  for(i = 0; i < len; ++i)
    err += res[i]*res[i];
  //err = sqrt(err) / len;
  if(err < tol)
  {
    M = 0;
    printf("No need to redistribute x with initial error=%lf.\n", err);
    number[10] = 0.0;
    number[11] = 0.0;
    charact[4] = 'O';
  }

  for(i = 0; i < len; ++i)
    p[i] = res[i];

  for(k = 1; k < M; ++k)
  {
    Ayfun(m, n, Ap, OMEGA, p);
    step_a = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < len; ++i)
    {
      step_a += res[i]*res[i];
      inner_prod += p[i]*Ap[i];
    }
    step_a = step_a/inner_prod;

    for(i = 1; i < n; ++i)
      for(j = 1; j < m; ++j)
	y[j][i] += step_a * p[i*m -i -m +j];

    for(i = 0; i < len; ++i)
      res_tilde[i] = res[i] - step_a * Ap[i];

    step_b = 0.0;
    inner_prod = 0.0;
    for(i = 0; i < len; ++i)
    {
      step_b += res_tilde[i]*res_tilde[i];
      inner_prod += res[i]*res[i];
    }
    step_b = step_b/inner_prod;

    for(i = 0; i < len; ++i)
    {
      res[i] = res_tilde[i];
      p[i] = res[i] + step_b * p[i];
    }

    err = 0.0;
    for(i = 0; i < len; ++i)
       err += res[i]*res[i];

    //err = sqrt(err) / len;
    if(err < tol)
      break;
  }

  if(M)
  {
    number[10] = (double)k;
    number[11] = err;
    if(err < tol)
    {
      printf("Mesh redistribution for x complete.\n");
      printf("\tNumber of ietrations: %d.\n\tTotal error :%g.\n", k, err);
      charact[4] = '$';
    }
    else
    {
      printf("Maximum number of iteration for x reached!\nThe total error is still %g\n", err);
      charact[4] = '#';
    }
  }


  free(res);
  free(res_tilde);
  free(p);
  free(Ap);
  free(xv);

  return 1;
}
