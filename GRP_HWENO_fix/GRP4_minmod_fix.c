#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>



#include "file_io.h"
#include "Riemann_solver.h"

#include "file_io_local.h"
#include "reconstruction.h"




int GRP4_minmod_fix
(double const CONFIG[], int const m, double const h,
 double rho[], double u[], double p[], runList *runhist, char *scheme)
{
  delete_runHist(runhist);
  int state;
  int i = 0, j = 0, k = 1, it = 0;

  char scheme_local[50] = "G4m2\0";
  printf("===========================\n");
  printf("The scheme [%s] started.\n", scheme_local);
  int len = 0;
  while(scheme_local[len] != '\0')
    ++len;
  ++len;
  for(k = 0; k < len; ++k)
    scheme[k] = scheme_local[k];


  clock_t tic, toc;
  int num_print_1 = 0, num_print_2 = 0, it_print;
  double current_cpu_time = 0.0, sum_cpu_time = 0.0, T = 0.0;

  
  
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
  double const thickness = CONFIG[10];
  int const    MaxStp    = (int)(CONFIG[11]);  // the number of time steps
  double const TIME      = CONFIG[12];
  int const    bod       = (int)CONFIG[14];
  int const    Primative = (int)CONFIG[15];
  int const    Deri      = (int)CONFIG[16];
  int const    Limiter   = (int)CONFIG[17];
  int const    Decomp    = (int)CONFIG[18];

  double running_info[N_RUNNING];
  running_info[3] = CONFIG[14];    // the boundary condition
  running_info[4] = CONFIG[16];    // the choice of the direvative reconstruction
  running_info[5] = CONFIG[18];    // use the charactoristic decomposition or not
  running_info[6] = CONFIG[17];    // use the limiter or not
  running_info[7] = CONFIG[9]; // threshold


  double c, stmp, D[4], U[4], wave_speed[2];


  double mom[m], ene[m];
  double rho_star, u_star, p_star, E_star;
  double rhoI[m+1], momI[m+1], eneI[m+1], uI[m+1], pI[m+1], half_rhoI[m+1], half_momI[m+1], half_eneI[m+1], half_uI[m+1], half_pI[m+1];
  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];

  double half_rho[m], half_u[m], half_p[m], half_mom[m], half_ene[m];

  double f01[m+1], f02[m+1], f03[m+1];
  double f11[m+1], f12[m+1], f13[m+1];
  double g01[m+1], g02[m+1], g03[m+1];
  double g11[m+1], g12[m+1], g13[m+1];
  double F1[m+1], F2[m+1], F3[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
			     * speed at each time step
			     */
  double tau, half_tau, nu, half_nu, alp, bet;


  double D0 = 1.0/6.0, D1 = 0.5-D0;


  if(Primative)
    for(j = 0; j < m; ++j)
    {
      mom[j] = rho[j]*u[j];
      ene[j] = p[j]/(gamma-1.0)+0.5*mom[j]*u[j];
    }
  else
    for(j = 0; j < m; ++j)
    {
      mom[j] = u[j];
      ene[j] = p[j];
      u[j] = mom[j]/rho[j];
      p[j] = (ene[j]-0.5*mom[j]*u[j])*(gamma-1.0);
    }

  running_info[0] = 0.0;
  running_info[1] = 0.0;
  running_info[2] = 0.0;
  WENO_50(running_info, m, h, eps, alp2, gamma, rho, mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);
//------------THE MAIN LOOP-------------
  for(k = 1; k <= MaxStp; ++k)
  {
    if(T+eps > TIME)
      break;


    state = insert_runHist(runhist);
    if(state)
    {
      printf("Not enough memory for the runhist node!\n\n");
      exit(100);//remains to modify the error code.
    }
    state = locate_runHist(k-1, runhist);
    if(state)
    {
      printf("The record has only %d compunonts while trying to reach runhist[%d].\n\n", state-1, k);
      exit(100);//remains to modify the error code.
    }


    running_info[0] = k;
    
    speed_max = 0.0;
    for(j = 0; j < m; ++j)
    {
      c = sqrt(gamma * p[j] / rho[j]);
      sigma = fabs(c) + fabs(u[j]);
      if(speed_max < sigma)
	speed_max = sigma;
    }
    tau = (CFL * h) / speed_max;
    if(T+tau > TIME){tau = TIME-T; T = TIME;} else{T += tau;}
    half_tau = 0.5*tau;
    nu = tau/h;
    half_nu = 0.5*nu;
    runhist->current->time[0] = tau;


    for(it_print = 0; it_print < num_print_1; ++it_print)
      printf(" ");
    printf("\r");
    fflush(stdout);
    num_print_1 = printf("%d | %g : %g | %g : %g", k, tau, T, current_cpu_time, sum_cpu_time);
    fflush(stdout);
    printf("\r");
    tic = clock();


    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			D_rho_L[j], D_u_L[j], 0.0, D_p_L[j],
			D_rho_R[j], D_u_R[j], 0.0, D_p_R[j],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      f01[j] = U[0]*U[1];
      f02[j] = f01[j]*U[1] + U[3];
      f03[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f01[j]*U[1]) * U[1];

      g01[j] = U[0]*D[1] + U[1]*D[0];
      g02[j] = D[0]*U[1]*U[1] + 2.0*U[0]*U[1]*D[1] + D[3];
      g03[j] = (D[3]*U[1] + U[3]*D[1])*gamma/(gamma-1.0);
      g03[j] = g03[j] + 0.5*D[0]*U[1]*U[1]*U[1] + 1.5*U[0]*U[1]*U[1]*D[1];

      rhoI[j] = U[0];
      momI[j] = U[0]*U[1];
      eneI[j] = 0.5*U[0]*U[1]*U[1] + U[3]/(gamma-1.0);
      half_rhoI[j] = rhoI[j] + half_tau*D[0];
      half_momI[j] = momI[j] + half_tau*(D[1]*U[0] + D[0]*U[1]);
      half_eneI[j] = eneI[j] + half_tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      uI[j] = U[1];
      pI[j] = U[3];
      half_uI[j] = U[1] + half_tau*D[1];
      half_pI[j] = U[3] + half_tau*D[3];
    }
    for(j = 0; j < m; ++j)
    {
      half_rho[j] = rho[j] - half_nu*((f01[j+1]-f01[j]) + 0.5*half_tau*(g01[j+1]-g01[j]));
      half_mom[j] = mom[j] - half_nu*((f02[j+1]-f02[j]) + 0.5*half_tau*(g02[j+1]-g02[j]));
      half_ene[j] = ene[j] - half_nu*((f03[j+1]-f03[j]) + 0.5*half_tau*(g03[j+1]-g03[j]));

      half_u[j] = half_mom[j] / half_rho[j];
      half_p[j] = (half_ene[j] - 0.5*half_mom[j]*half_u[j]) * (gamma-1.0);
      if(half_p[j] < 0.0)
        printf("half(%d,%d), %g\n", k, j, half_p[j]);
    }

    running_info[1] = T - half_tau;  // time
    running_info[2] = 1.0;           // half
    GRP_minmod(running_info, m, h, alp2, half_rho, half_u, half_p, half_rhoI, half_uI, half_pI, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			D_rho_L[j], D_u_L[j], 0.0, D_p_L[j],
			D_rho_R[j], D_u_R[j], 0.0, D_p_R[j],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      g11[j] = U[0]*D[1] + U[1]*D[0];
      g12[j] = D[0]*U[1]*U[1] + 2.0*U[0]*U[1]*D[1] + D[3];
      g13[j] = (D[3]*U[1] + U[3]*D[1])*gamma/(gamma-1.0);
      g13[j] = g13[j] + 0.5*D[0]*U[1]*U[1]*U[1] + 1.5*U[0]*U[1]*U[1]*D[1];

      F1[j] = f01[j] + tau*(D0*g01[j] + D1*g11[j]);
      F2[j] = f02[j] + tau*(D0*g02[j] + D1*g12[j]);
      F3[j] = f03[j] + tau*(D0*g03[j] + D1*g13[j]);

      rhoI[j] = rhoI[j] + tau*D[0];
      momI[j] = momI[j] + tau*(D[1]*U[0] + D[0]*U[1]);
      eneI[j] = eneI[j] + tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      uI[j] = uI[j] + tau*D[1];
      pI[j] = pI[j] + tau*D[3];
    }

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
	rho[j] = rho[j] - nu*(F1[j+1]-F1[j]);
	mom[j] = mom[j] - nu*(F2[j+1]-F2[j]);
	ene[j] = ene[j] - nu*(F3[j+1]-F3[j]);
	u[j] = mom[j] / rho[j];
	p[j] = (ene[j] - 0.5*mom[j]*u[j])*(gamma-1.0);

      if(p[j] < 0.0)
        printf("    (%d,%d), %g\n", k, j, p[j]);
    }

    running_info[1] = T;
    running_info[2] = 0.0;
    GRP_minmod(running_info, m, h, alp2, rho, u, p, rhoI, uI, pI, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);

    toc = clock();
    runhist->current->time[1] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;
    current_cpu_time = runhist->current->time[1];
    sum_cpu_time += runhist->current->time[1];
  }
  k = k-1;

  if(check_runHist(runhist))
  {
    printf("The runhist->length is %d.\nBut the number of the runNodes are %d.\n\n", runhist->length, runhist->length - check_runHist(runhist));
    exit(100);
  }
  if(k - runhist->length)
  {
    printf("After %d steps of computation, the number of the runNodes are %d.\n\n", k, runhist->length);
    exit(100);
  }


  printf("The cost of CPU time for [%s] computing this\n", scheme);
  printf("problem to time %g with %d steps is %g seconds.\n", T, k, sum_cpu_time);
  printf("===========================\n");
//------------END OFQ THE MAIN LOOP-------------


  return k;
}
