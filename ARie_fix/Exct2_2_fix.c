#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "Riemann_solver.h"

#ifndef L_STR
#include "file_io_local.h"
#endif

#include "reconstruction.h"




int Exct2_2_fix
(double const CONFIG[], int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label)
{
  delete_runHist(runhist);
  int i = 0, j = 0, k = 1, it = 0;
  int state, len = 0;
  char scheme[L_STR] = "Exct2-2\0";
  char version[L_STR], err_msg[L_STR];
  strcpy(version, add_mkdir);
  strcpy(add_mkdir, "../SOLUTION/\0");
  state = make_directory(add_mkdir, err_msg, label, scheme, version, m, 1, CONFIG);
  if(state)
  {
    printf("%s", err_msg);
    exit(state);
  }


  printf("===========================\n");
  printf("The scheme [%s] started.\n", scheme);


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
  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];
  double rho1[m], mom1[m], ene1[m], u1[m], p1[m];
  double uL, uR, pL, pR, DuL, DuR, DpL, DpR;

  double F1[m+1], F2[m+1], F3[m+1], rhoI[m+1], uI[m+1], pI[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
					 * speed at each time step
					 */
  double tau, alp, bet, nu;


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
  running_info[1] = 0.0;
  GRP_minmod0(running_info, m, h, alp2, rho, mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);
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
    
    
    running_info[0] = (double)k;

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
    nu = tau/h;
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
      uL = u_L[j]/rho_L[j];                           uR = u_R[j]/rho_R[j];
      pL = (gamma-1.0)*(p_L[j] - 0.5*rho_L[j]*uL*uL); pR = (gamma-1.0)*(p_R[j] - 0.5*rho_R[j]*uR*uR);
      DuL = (D_u_L[j] - D_rho_L[j]*uL)/rho_L[j];      DuR = (D_u_R[j] - D_rho_R[j]*uR)/rho_R[j];
      DpL = (gamma-1.0)*(D_p_L[j] - 0.5*D_rho_L[j]*uL*uL - u_L[j]*DuL);
      DpR = (gamma-1.0)*(D_p_R[j] - 0.5*D_rho_R[j]*uR*uR - u_R[j]*DuR);
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], uL, 0.0, pL,
			rho_R[j], uR, 0.0, pR,
			D_rho_L[j], DuL, 0.0, DpL,
			D_rho_R[j], DuR, 0.0, DpR,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      F1[j] = U[0]*U[1];
      F2[j] = U[0]*U[1]*U[1]  + U[3];
      F3[j] = U[3]*U[1]*gamma/(gamma-1.0);
      F3[j] = 0.5*U[0]*U[1]*U[1]*U[1] + F3[j];
    }

    for(j = 0; j < m; ++j)
    {
      rho1[j] = rho[j] - nu*(F1[j+1]-F1[j]);
      mom1[j] = mom[j] - nu*(F2[j+1]-F2[j]);
      ene1[j] = ene[j] - nu*(F3[j+1]-F3[j]);

	u1[j] = mom1[j] / rho1[j];
	p1[j] = (ene1[j] - 0.5*mom1[j]*u1[j])*(gamma-1.0);
    }

    running_info[1] = T;
    running_info[2] = 1.0;
    GRP_minmod0(running_info, m, h, alp2, rho1, mom1, ene1, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


    for(j = 0; j < m+1; ++j)
    {
      uL = u_L[j]/rho_L[j];                           uR = u_R[j]/rho_R[j];
      pL = (gamma-1.0)*(p_L[j] - 0.5*rho_L[j]*uL*uL); pR = (gamma-1.0)*(p_R[j] - 0.5*rho_R[j]*uR*uR);
      DuL = (D_u_L[j] - D_rho_L[j]*uL)/rho_L[j];      DuR = (D_u_R[j] - D_rho_R[j]*uR)/rho_R[j];
      DpL = (gamma-1.0)*(D_p_L[j] - 0.5*D_rho_L[j]*uL*uL - u_L[j]*DuL);
      DpR = (gamma-1.0)*(D_p_R[j] - 0.5*D_rho_R[j]*uR*uR - u_R[j]*DuR);
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], uL, 0.0, pL,
			rho_R[j], uR, 0.0, pR,
			D_rho_L[j], DuL, 0.0, DpL,
			D_rho_R[j], DuR, 0.0, DpR,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      F1[j] = U[0]*U[1];
      F2[j] = U[0]*U[1]*U[1]  + U[3];
      F3[j] = U[3]*U[1]*gamma/(gamma-1.0);
      F3[j] = 0.5*U[0]*U[1]*U[1]*U[1] + F3[j];
    }

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
      rho[j] = 0.5*((rho[j] + rho1[j]) - nu*(F1[j+1]-F1[j]));
      mom[j] = 0.5*((mom[j] + mom1[j]) - nu*(F2[j+1]-F2[j]));
      ene[j] = 0.5*((ene[j] + ene1[j]) - nu*(F3[j+1]-F3[j]));

	u[j] = mom[j] / rho[j];
	p[j] = (ene[j] - 0.5*mom[j]*u[j])*(gamma-1.0);
    }

    running_info[1] = T;
    running_info[2] = 0.0;
    GRP_minmod0(running_info, m, h, alp2, rho, mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


    toc = clock();
    runhist->current->time[1] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;
    current_cpu_time = runhist->current->time[1];
    sum_cpu_time += runhist->current->time[1];
  }
  k = k-1;
  if(check_runHist(runhist))
  {
    printf("The runhist->length is %d.\nBut the number of the runNodes is %d.\n\n", runhist->length, runhist->length - check_runHist(runhist));
    exit(100);
  }
  if(k - runhist->length)
  {
    printf("After %d steps of computation, the number of the runNodes is %d.\n\n", k, runhist->length);
    exit(100);
  }
  

  printf("The cost of CPU time for [%s] computing this\n", scheme);
  printf("problem to time %g with %d steps is %g seconds.\n", T, k, sum_cpu_time);
  printf("===========================\n");
//------------END OFQ THE MAIN LOOP-------------


  return k;
}
