#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "file_io_2D.h"

#include "solver.h"
#include "reconstruction.h"





/* This function use Godunov scheme to solve 1-D
 * scalar conservation law:
 *                    u_t + f(u)_x = 0.
 *U      is the gird function aproximating the solution u.
 *time   denotes the value of each t^n since tau is
 *         changable.
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *m      is the number of the grids.
 */
int FV_WENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList *runhist, char *scheme)
{
  int i = 0, j = 0, k = 1, it = 0;  /* j is a frequently used index for
				     * spatial variables. n is a frequ-
				     * ently used index for the time
				     * step.
				     */
  char scheme_local[50] = "WENO-5-RF-4\0";
  printf("===========================\n");
  printf("The scheme [%s] started.\n", scheme_local);
  int len = 0;
  while(scheme_local[len] != '\0')
    ++len;
  ++len;
  for(k = 0; k < len; ++k)
    scheme[k] = scheme_local[k];


  clock_t tic, toc;
  double sum = 0.0, T = 0.0;
  int const    MaxStp     = (int)(OPT[0]);  // the number of time steps
  double const TIME       = OPT[1];
  int const    inter_data = (int)OPT[2];
  int const    bod        = (int)OPT[4];
  int const    Riemann    = (int)OPT[5];
  int const    WENOD      = (int)OPT[6];
  int const    limiter    = (int)OPT[7];
  int vk0, vk1;

  int running_info[5];
  running_info[2] = bod;
  running_info[3] = WENOD;
  running_info[4] = limiter;

  double const gamma = CONFIG[0];  // the constant of the perfect gas
  double const CFL = CONFIG[1];    // CFL number
  double const eps = CONFIG[2];    // the largest value could be treat as zero
  double const alp1 = CONFIG[3];
  double const alp2 = CONFIG[4];
  double const bet1 = CONFIG[5];
  double const bet2 = CONFIG[6];
  double const modifier = CONFIG[7];
  double const tol = CONFIG[8];

  double c, stmp, D[4], U[4], wave_speed[2];

  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];


  double mom[m], ene[m];

  double rho_1[m], mom_1[m], ene_1[m];
  double rho_2[m], mom_2[m], ene_2[m];
  double rho_3[m], mom_3[m], ene_3[m];


  double f01[m+1], f02[m+1], f03[m+1];
  double f11[m+1], f12[m+1], f13[m+1];
  double f21[m+1], f22[m+1], f23[m+1];
  double f31[m+1], f32[m+1], f33[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
					 * speed at each time step
					 */
  double tau, half_tau, alp, bet;


double a10 = 1.0;
double a20 = 649.0/1600.0, a21 = 951.0/1600.0;
double a30 = 53989.0/2500000.0, a31 = 4806213.0/20000000.0, a32 = 23619.0/32000.0;
double a40 = 0.2, a41 = 6127.0/30000.0, a42 = 7873.0/30000.0, a43=1.0/3.0;
double b10 = 0.5;
double b20 = -10890423.0/25193600.0, b21 = 5000.0/7873.0;
double b30 = -102261.0/5000000.0, b31 = -5121.0/20000.0, b32 = 7873.0/10000.0;
double b40 = 0.1, b41 = 1.0/6.0, b42 = 0.0, b43 = 1.0/6.0;



  for(j = 0; j < m; ++j)
  {
    mom[j] = rho[0][j]*u[0][j];
    ene[j] = p[0][j]/(gamma-1.0)+0.5*mom[j]*u[0][j];
  }

//------------THE MAIN LOOP-------------
  for(k = 1; k <= MaxStp; ++k)
  {    
    if(T+eps > TIME)
      break;


    insert_runList(runhist);
    if(!runhist->tail)
    {
      printf("Not enough memory for the runhist node!\n\n");
      exit(100);
    }
    locate_runList(k, runhist);

    vk0 = (k-1)*inter_data;  // vk0==0 or vk0==k-1
    vk1 =     k*inter_data;  // vk1==0 or vk1==k
    running_info[0] = k;
    //printf("-----------------%d-----------------", k);
    speed_max = 0.0;
    for(j = 0; j < m; ++j)
      {
	c = sqrt(gamma * p[vk0][j] / rho[vk0][j]);
	sigma = fabs(c) + fabs(u[vk0][j]);
	speed_max = ((speed_max < sigma) ? sigma : speed_max);
      }
    tau = (CFL * h) / speed_max;
    if(T+tau > TIME){tau = TIME-T; T = TIME;} else{T += tau;}
    half_tau = 0.5*tau;
    runhist->current->time[0] = tau;
    //printf("%g, %g\n", tau, T);

    tic = clock();


    //======FIRST=====================
    running_info[1] = 0;

    WENO_5_noD(running_info, m, h, eps, alp2, gamma, rho[vk0], mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R);
    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, rho_L[j], rho_R[j], 0.0, 0.0, u_L[j], u_R[j], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_L[j], p_R[j], 0.0, 0.0, gamma, eps);

      f01[j] = U[0]*U[1];
      f02[j] = f01[j]*U[1] + U[3];
      f03[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f01[j]*U[1]) * U[1];
    }

    for(j = 0; j < m; ++j)
    {
      rho_1[j] = rho[vk0][j] - half_tau*(f01[j+1]-f01[j])/h;
      mom_1[j] =      mom[j] - half_tau*(f02[j+1]-f02[j])/h;
      ene_1[j] =      ene[j] - half_tau*(f03[j+1]-f03[j])/h;
    }

    //======SECOND=====================
    running_info[1] = 1;

    WENO_5_noD(running_info, m, h, eps, alp2, gamma, rho_1, mom_1, ene_1, rho_L, rho_R, u_L, u_R, p_L, p_R);
    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, rho_L[j], rho_R[j], 0.0, 0.0, u_L[j], u_R[j], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_L[j], p_R[j], 0.0, 0.0, gamma, eps);

      f11[j] = U[0]*U[1];
      f12[j] = f11[j]*U[1] + U[3];
      f13[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f11[j]*U[1]) * U[1];
    }

    for(j = 0; j < m; ++j)
    {
      rho_2[j] = (a20*rho[vk0][j]+a21*rho_1[j]) - tau*(b20*(f01[j+1]-f01[j])+b21*(f11[j+1]-f11[j]))/h;
      mom_2[j] = (a20*     mom[j]+a21*mom_1[j]) - tau*(b20*(f02[j+1]-f02[j])+b21*(f12[j+1]-f12[j]))/h;
      ene_2[j] = (a20*     ene[j]+a21*ene_1[j]) - tau*(b20*(f03[j+1]-f03[j])+b21*(f13[j+1]-f13[j]))/h;
    }

    //======THIRD=====================
    running_info[1] = 2;

    WENO_5_noD(running_info, m, h, eps, alp2, gamma, rho_2, mom_2, ene_2, rho_L, rho_R, u_L, u_R, p_L, p_R);
    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, rho_L[j], rho_R[j], 0.0, 0.0, u_L[j], u_R[j], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_L[j], p_R[j], 0.0, 0.0, gamma, eps);

      f21[j] = U[0]*U[1];
      f22[j] = f21[j]*U[1] + U[3];
      f23[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f21[j]*U[1]) * U[1];
    }

    for(j = 0; j < m; ++j)
    {
      rho_3[j] = (a30*rho[vk0][j]+a31*rho_1[j]+a32*rho_2[j]) - tau*(b30*(f01[j+1]-f01[j])+b31*(f11[j+1]-f11[j])+b32*(f21[j+1]-f21[j]))/h;
      mom_3[j] = (a30*     mom[j]+a31*mom_1[j]+a32*mom_2[j]) - tau*(b30*(f02[j+1]-f02[j])+b31*(f12[j+1]-f12[j])+b32*(f22[j+1]-f22[j]))/h;
      ene_3[j] = (a30*     ene[j]+a31*ene_1[j]+a32*ene_2[j]) - tau*(b30*(f03[j+1]-f03[j])+b31*(f13[j+1]-f13[j])+b32*(f23[j+1]-f23[j]))/h;
    }

    //======FORTH=====================
    running_info[1] = 3;

    WENO_5_noD(running_info, m, h, eps, alp2, gamma, rho_3, mom_3, ene_3, rho_L, rho_R, u_L, u_R, p_L, p_R);
    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, rho_L[j], rho_R[j], 0.0, 0.0, u_L[j], u_R[j], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, p_L[j], p_R[j], 0.0, 0.0, gamma, eps);

      f31[j] = U[0]*U[1];
      f32[j] = f31[j]*U[1] + U[3];
      f33[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f31[j]*U[1]) * U[1];
    }




//===============THE CORE ITERATION=================
//*
    for(j = 0; j < m; ++j)
    {
      rho[vk1][j] = (a40*rho[vk0][j]+a41*rho_1[j]+a42*rho_2[j]+a43*rho_3[j]) - tau*(b40*(f01[j+1]-f01[j])+b41*(f11[j+1]-f11[j])+b42*(f21[j+1]-f21[j])+b43*(f31[j+1]-f31[j]))/h;
           mom[j] = (a40*     mom[j]+a41*mom_1[j]+a42*mom_2[j]+a43*mom_3[j]) - tau*(b40*(f02[j+1]-f02[j])+b41*(f12[j+1]-f12[j])+b42*(f22[j+1]-f22[j])+b43*(f32[j+1]-f32[j]))/h;
           ene[j] = (a40*     ene[j]+a41*ene_1[j]+a42*ene_2[j]+a43*ene_3[j]) - tau*(b40*(f03[j+1]-f03[j])+b41*(f13[j+1]-f13[j])+b42*(f23[j+1]-f23[j])+b43*(f33[j+1]-f33[j]))/h;

      u[vk1][j] = mom[j] / rho[vk1][j];
      p[vk1][j] = (ene[j] - 0.5*mom[j]*u[vk1][j])*(gamma-1.0);
    }
    /*/
    for(j = 0; j < m; ++j)
    {
      rho[vk1][j] = rho[vk0][j] - tau*(f01[j+1]-f01[j])/h;
           mom[j] =      mom[j] - tau*(f02[j+1]-f02[j])/h;
           ene[j] =      ene[j] - tau*(f03[j+1]-f03[j])/h;

      u[vk1][j] = mom[j] / rho[vk1][j];
      p[vk1][j] = (ene[j] - 0.5*mom[j]*u[vk1][j])*(gamma-1.0);
    }
    //*/

    toc = clock();
    runhist->current->time[1] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += runhist->current->time[1];
  }
  k = k-1;

  if(check_runList(runhist))
  {
    printf("The runhist->length is %d.\nBut the number of the runNodes are %d.\n\n", runhist->length, runhist->length - check_runList(runhist));
    exit(100);
  }
  if(k - runhist->length)
  {
    printf("After %d steps of computation, the number of the runNodes are %d.\n\n", k, runhist->length);
    exit(100);
  }


  printf("The cost of CPU time for [%s] solving this problem by %d steps is %g seconds.\n", scheme, k, sum);
  printf("===========================\n");
//------------END OFQ THE MAIN LOOP-------------


  return k;
}
