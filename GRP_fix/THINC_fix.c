/**************************************************************************
 * This function use the GRP scheme to solve 1-D Euler eqautions:
 *                    w_t + f(w)_x = 0
 * on fixed meshes, where
 *      / rho  \       /   rho*u   \
 * w = | rho*u |  f = | rho*u^2+p  |   E=0.5*rho*u^2 + p/(gamma-1)
 *     \   E  / ,     \ rho*(E+p) / ,                             .
 *
 *
 * CONFIG[0] is the constant of the perfect gas
 * CONFIG[1] is the CFL number
 * CONFIG[2] is the largest value can be seen as zero
 * CONFIG[3] is the first limiter of the slope
 * CONFIG[4] is the second limiter of the slope
 * CONFIG[5] is the first parameter of the monitor function
 * CONFIG[6] is the second parameter of the monitor function
 * CONFIG[7] is the modifier of the mesh redistribution
 * CONFIG[8] is the tolerance of the mesh redistribution
 *
 * OPT[0] is the maximal step to compute.
 * OPT[1] is the time to stop the computation
 * OPT[2] is the switch of whether keep the inter-data during the computation
 * OPT[3] is the switch of wether use an adaptive mesh
 * OPT[4] controls the choice of the boundary condition
 *
 * m    is the number of the spatial grids
 * h    is the size of the spatial grids
 * rho  is the density
 * u    is the velocity
 * p    is the pressure
 * runhist_n[]
 * runhist_c
 **************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "file_io.h"
#include "Riemann_solver.h"

#include "file_io_local.h"
#include "reconstruction.h"




int THINC_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList *runhist, char *scheme)
{
  int i = 0, j = 0, k = 1, it = 0;  /* j is a frequently used index for
				     * spatial variables. n is a frequ-
				     * ently used index for the time
				     * step.
				     */
  char scheme_local[50] = "G2m2\0";
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
  int vk0, vk1;
  
  
  double const gamma     = CONFIG[0];  // the constant of the perfect gas
  double const CFL       = CONFIG[1];  // CFL number
  double const eps       = CONFIG[2];  // the largest value could be treat as zero
  double const alp1      = CONFIG[3];
  double const alp2      = CONFIG[4];
  double const bet1      = CONFIG[5];
  double const bet2      = CONFIG[6];
  double const modifier  = CONFIG[7];
  double const tol       = CONFIG[8];
  double const threshold = CONFIG[9];

  int const    MaxStp     = (int)(OPT[0]);  // the number of time steps
  double const TIME       = OPT[1];
  int const    inter_data = (int)OPT[2];
  int const    bod        = (int)OPT[4];
  int const    Riemann    = (int)OPT[5];
  int const    WENOD      = (int)OPT[6];
  int const    decomp     = (int)OPT[7];
  int const    limiter    = (int)OPT[8];

  double running_info[N_RUNNING];
  running_info[3] = OPT[4];    // the boundary condition
  running_info[4] = OPT[6];    // the choice of the direvative reconstruction
  running_info[5] = OPT[7];    // use the charactoristic decomposition or not
  running_info[6] = OPT[8];    // use the limiter or not
  running_info[7] = CONFIG[9]; // threshold
  int trouble[m];
  for(j = 0; j < m; ++j)
    trouble[j] = 0;


  double c, stmp, D[4], U[4], wave_speed[2];


  double mom[m], ene[m];
  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];

  double F1[m+1], F2[m+1], F3[m+1], rhoI[m+1], uI[m+1], pI[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
					 * speed at each time step
					 */
  double tau, half_tau, alp, bet, nu;


  if(Riemann)
    for(j = 0; j < m; ++j)
    {
      mom[j] = rho[0][j]*u[0][j];
      ene[j] = p[0][j]/(gamma-1.0)+0.5*mom[j]*u[0][j];
    }
  else
    for(j = 0; j < m; ++j)
    {
      mom[j] = u[0][j];
      ene[j] = p[0][j];
      u[0][j] = mom[j]/rho[0][j];
      p[0][j] = (ene[j]-0.5*mom[j]*u[0][j])*(gamma-1.0);
    }


  running_info[0] = 0.0;
  running_info[1] = 0.0;
  running_info[1] = 0.0;
  THINC0(running_info, m, h, alp2, rho[0], u[0], p[0], rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R, trouble);

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
    running_info[0] = (double)k;
    //printf("-----------------%d-----------------", k);
    speed_max = 0.0;
    for(j = 0; j < m; ++j)
    {
      c = sqrt(gamma * p[vk0][j] / rho[vk0][j]);
      sigma = fabs(c) + fabs(u[vk0][j]);
      if(speed_max < sigma)
        speed_max = sigma;
    }
    tau = (CFL * h) / speed_max;

    if(T+tau > TIME){tau = TIME-T; T = TIME;} else{T += tau;}
    half_tau = 0.5*tau;
    nu = tau/h;
    runhist->current->time[0] = tau;

    tic = clock();


    for(j = 0; j < m+1; ++j)
      {/*
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			D_rho_L[j], D_u_L[j], 0.0, D_p_L[j],
			D_rho_R[j], D_u_R[j], 0.0, D_p_R[j],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);/*/
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);//*/

      F1[j] = U[0]*U[1] + half_tau*(D[0]*U[1]+U[0]*D[1]);
      F2[j] = U[0]*U[1]*U[1] + half_tau*(D[0]*U[1]*U[1]+2.0*U[0]*U[1]*D[1]) + (U[3]+half_tau*D[3]);
      F3[j] = (U[3]*U[1]+half_tau*(D[3]*U[1]+U[3]*D[1]))*gamma/(gamma-1.0);
      F3[j] = 0.5*U[0]*U[1]*U[1]*U[1] + half_tau*(0.5*D[0]*U[1]*U[1]*U[1]+1.5*U[0]*U[1]*U[1]*D[1]) + F3[j];

      rhoI[j] = U[0] + tau*D[0];
        uI[j] = U[1] + tau*D[1];
        pI[j] = U[3] + tau*D[3];
    }

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
      rho[vk1][j] = rho[vk0][j] - nu*(F1[j+1]-F1[j]);
	   mom[j] =      mom[j] - nu*(F2[j+1]-F2[j]);
	   ene[j] =      ene[j] - nu*(F3[j+1]-F3[j]);

	u[vk1][j] = mom[j] / rho[vk1][j];
	p[vk1][j] = (ene[j] - 0.5*mom[j]*u[vk1][j])*(gamma-1.0);
    } 
//==================================================
    write_column(m, rho[vk1], "rho_sheer", "running");
    write_column(m, u[vk1], "u_sheer", "running");
    write_column(m, p[vk1], "p_sheer", "running");
    write_column(m, mom, "momx_sheer", "running");
    write_column(m, ene, "ene_sheer", "running");

    running_info[1] = T;
    running_info[2] = 0.0;
    THINC0(running_info, m, h, alp2, rho[vk1], u[vk1], p[vk1], rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R, trouble);


    toc = clock();
    runhist->current->time[1] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += runhist->current->time[1];
  }
  k = k-1;
  if(check_runList(runhist))
  {
    printf("The runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runList(runhist));
    exit(100);
  }
  if(k - runhist->length)
  {
    printf("After %d steps of computation, the number of the runNodes is %d.\n\n", k, runhist->length);
    exit(100);
  }
  

  printf("The cost of CPU time for [%s] solving this problem by %d steps is %g seconds.\n", scheme, k, sum);
  printf("===========================\n");
//------------END OFQ THE MAIN LOOP-------------


  return k;
}
