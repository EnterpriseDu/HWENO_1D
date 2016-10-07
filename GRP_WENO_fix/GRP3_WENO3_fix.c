#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>



#include "file_io.h"
#include "Riemann_solver.h"

#include "file_io_local.h"
#include "reconstruction.h"


int GRP3_WENO3_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList *runhist, char *scheme)
{
  int i = 0, j = 0, k = 1, it = 0;

  char scheme_local[50] = "G3W3\0";
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
  double const thickness = CONFIG[10];

  int const    MaxStp     = (int)(OPT[0]);  // the number of time steps
  double const TIME       = OPT[1];
  int const    inter_data = (int)OPT[2];
  int const    bod        = (int)OPT[4];
  int const    Riemann    = (int)OPT[5];
  int const    WENOD      = (int)OPT[6];
  int const    decomp     = (int)OPT[7];
  int const    limiter    = (int)OPT[8];

  double running_info[N_RUNNING];
  running_info[3] = OPT[4];  // the boundary condition
  running_info[4] = OPT[6];  // the choice of the direvative reconstruction
  running_info[5] = OPT[7];  // use the charactoristic decomposition or not
  running_info[6] = OPT[8];    // use the limiter or not
  running_info[7] = CONFIG[9]; // threshold
  running_info[8] = CONFIG[10]; //thickness



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
  double tau, nu, alp, bet;




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

  running_info[0] = 0.0;  // k
  running_info[1] = 0.0;  // time
  running_info[2] = 0.0;  // not half
  WENO_30(running_info, m, h, eps, alp2, gamma, rho[0], mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);
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
    runhist->current->trouble0 = (int *)malloc(sizeof(int) * m);
    runhist->current->trouble1 = (int *)malloc(sizeof(int) * m);
    if((!runhist->current->trouble0) || (!runhist->current->trouble1))
    {
      printf("Not enough memory for the runhist node!\n\n");
      exit(100);
    }

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
    nu = tau/h;
    runhist->current->time[0] = tau;
    //printf("%g, %g\n", tau, T);

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

      //rhoI[j] = U[0];
      //momI[j] = U[0]*U[1];
      //eneI[j] = 0.5*U[0]*U[1]*U[1] + U[3]/(gamma-1.0);
      //half_rhoI[j] = rhoI[j] + tau*D[0];
      //half_momI[j] = momI[j] + tau*(D[1]*U[0] + D[0]*U[1]);
      //half_eneI[j] = eneI[j] + tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      //uI[j] = U[1];
      //pI[j] = U[3];
      //half_uI[j] = U[1] + tau*D[1];
      //half_pI[j] = U[3] + tau*D[3];
    }
    for(j = 0; j < m; ++j)
    {
      half_rho[j] = rho[vk0][j] - nu*((f01[j+1]-f01[j]) + 0.5*tau*(g01[j+1]-g01[j]));
      half_mom[j] =      mom[j] - nu*((f02[j+1]-f02[j]) + 0.5*tau*(g02[j+1]-g02[j]));
      half_ene[j] =      ene[j] - nu*((f03[j+1]-f03[j]) + 0.5*tau*(g03[j+1]-g03[j]));

      half_u[j] = half_mom[j] / half_rho[j];
      half_p[j] = (half_ene[j] - 0.5*half_mom[j]*half_u[j]) * (gamma-1.0);
      if(half_p[j] < 0.0)
        printf("half(%d,%d), %g\n", k, j, half_p[j]);
    }

    running_info[1] = T;  // time
    running_info[2] = 1.0;           // half
    WENO_30(running_info, m, h, eps, alp2, gamma, half_rho, half_mom, half_ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			D_rho_L[j], D_u_L[j], 0.0, D_p_L[j],
			D_rho_R[j], D_u_R[j], 0.0, D_p_R[j],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      f11[j] = U[0]*U[1];
      f12[j] = f11[j]*U[1] + U[3];
      f13[j] = (gamma*U[3]/(gamma-1.0) + 0.5*f11[j]*U[1]) * U[1];

      g11[j] = U[0]*D[1] + U[1]*D[0];
      g12[j] = D[0]*U[1]*U[1] + 2.0*U[0]*U[1]*D[1] + D[3];
      g13[j] = (D[3]*U[1] + U[3]*D[1])*gamma/(gamma-1.0);
      g13[j] = g13[j] + 0.5*D[0]*U[1]*U[1]*U[1] + 1.5*U[0]*U[1]*U[1]*D[1];

      F1[j] = (2.0*f01[j] + f11[j] + 0.5*tau*g11[j])/3.0;
      F2[j] = (2.0*f02[j] + f12[j] + 0.5*tau*g12[j])/3.0;
      F3[j] = (2.0*f03[j] + f13[j] + 0.5*tau*g13[j])/3.0;

      //rhoI[j] = rhoI[j] + tau*D[0];
      //momI[j] = momI[j] + tau*(D[1]*U[0] + D[0]*U[1]);
      //eneI[j] = eneI[j] + tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      //uI[j] = uI[j] + tau*D[1];
      //pI[j] = pI[j] + tau*D[3];
    }

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
      rho[vk1][j] = rho[vk0][j] - nu*(F1[j+1]-F1[j]);
           mom[j] =      mom[j] - nu*(F2[j+1]-F2[j]);
	   ene[j] =      ene[j] - nu*(F3[j+1]-F3[j]);
	u[vk1][j] = mom[j] / rho[vk1][j];
	p[vk1][j] = (ene[j] - 0.5*mom[j]*u[vk1][j])*(gamma-1.0);

      if(p[vk1][j] < 0.0)
        printf("    (%d,%d), %g\n", k, j, p[vk1][j]);
    }

    running_info[1] = T;
    running_info[2] = 0.0;  // not half
    WENO_30(running_info, m, h, eps, alp2, gamma, rho[vk1], mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);

    toc = clock();
    runhist->current->time[1] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;
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