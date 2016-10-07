#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>



#include "file_io.h"
#include "Riemann_solver.h"

#include "file_io_local.h"
#include "reconstruction.h"


int GRP5_WENO5_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList *runhist, char *scheme)
{
  int i = 0, j = 0, k = 1, it = 0;

  char scheme_local[50] = "G5W5\0";
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
  double rhoI[m+1], momI[m+1], eneI[m+1], uI[m+1], pI[m+1];//, half_rhoI[m+1], half_momI[m+1], half_eneI[m+1], half_uI[m+1], half_pI[m+1];
  double rhoI1[m+1], momI1[m+1], eneI1[m+1], uI1[m+1], pI1[m+1];
  double rhoI2[m+1], momI2[m+1], eneI2[m+1], uI2[m+1], pI2[m+1];
  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];

  double rho1[m], u1[m], p1[m], mom1[m], ene1[m];
  double rho2[m], u2[m], p2[m], mom2[m], ene2[m];

  double f01[m+1], f02[m+1], f03[m+1];
  double g01[m+1], g02[m+1], g03[m+1];
  double g11[m+1], g12[m+1], g13[m+1];
  double g21[m+1], g22[m+1], g23[m+1];
  double F11[m+1], F12[m+1], F13[m+1];
  double F21[m+1], F22[m+1], F23[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
			     * speed at each time step
			     */
  double tau, tau1, tau2, nu, nu1, nu2, alp, bet;


  double D10 = -0.5, D11 = 1.5;
  double D20 = 0.25, D21 = 25.0/36.0,  D22 = 1.0/18.0;


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
  WENO_50(running_info, m, h, eps, alp2, gamma, rho[0], mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);
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
    tau1 = 0.4*tau;
    nu = tau/h;
    nu1 = 0.4*nu;
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

      rhoI[j] = U[0];
      momI[j] = U[0]*U[1];
      eneI[j] = 0.5*U[0]*U[1]*U[1] + U[3]/(gamma-1.0);
      //rhoI1[j] = rhoI[j] + tau1*D[0];
      //momI1[j] = momI[j] + tau1*(D[1]*U[0] + D[0]*U[1]);
      //eneI1[j] = eneI[j] + tau1*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      uI[j] = U[1];
      pI[j] = U[3];
      //uI1[j] = U[1] + tau1*D[1];
      //pI1[j] = U[3] + tau1*D[3];
    }
    for(j = 0; j < m; ++j)
    {
      rho1[j] = rho[vk0][j] - nu1*((f01[j+1]-f01[j]) + 0.5*tau1*(g01[j+1]-g01[j]));
      mom1[j] =      mom[j] - nu1*((f02[j+1]-f02[j]) + 0.5*tau1*(g02[j+1]-g02[j]));
      ene1[j] =      ene[j] - nu1*((f03[j+1]-f03[j]) + 0.5*tau1*(g03[j+1]-g03[j]));

      u1[j] = mom1[j] / rho1[j];
      p1[j] = (ene1[j] - 0.5*mom1[j]*u1[j]) * (gamma-1.0);
      if(p1[j] < 0.0)
        printf("half(%d,%d), %g\n", k, j, p1[j]);
    }

    running_info[1] = T - tau + tau1;  // time
    running_info[2] = 1.0;           // half
    //WENO_5(running_info, m, h, eps, alp2, gamma, rho1, mom1, ene1, rhoI1, uI1, pI1, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R, runhist->current->trouble1);
    WENO_50(running_info, m, h, eps, alp2, gamma, rho1, mom1, ene1, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


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

      F11[j] = f01[j] + 0.5*tau*(D10*g01[j] + D11*g11[j]);
      F12[j] = f02[j] + 0.5*tau*(D10*g02[j] + D11*g12[j]);
      F13[j] = f03[j] + 0.5*tau*(D10*g03[j] + D11*g13[j]);

      //rhoI2[j] = rhoI[j] + tau*D[0];
      //momI2[j] = momI[j] + tau*(D[1]*U[0] + D[0]*U[1]);
      //eneI2[j] = eneI[j] + tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      //uI2[j] = uI[j] + tau*D[1];
      //pI2[j] = pI[j] + tau*D[3];
    }
    for(j = 0; j < m; ++j)
    {
      rho2[j] = rho[vk0][j] - nu*(F11[j+1]-F11[j]);
      mom2[j] =      mom[j] - nu*(F12[j+1]-F12[j]);
      ene2[j] =      ene[j] - nu*(F13[j+1]-F13[j]);
      //rho2[j] = rho[vk0][j] - nu*((f01[j+1]-f01[j]) + 0.5*tau*(D10*(g01[j+1]-g01[j])+D11*(g11[j+1]-g11[j])));
      //mom2[j] =      mom[j] - nu*((f02[j+1]-f02[j]) + 0.5*tau*(D10*(g02[j+1]-g02[j])+D11*(g12[j+1]-g12[j])));
      //ene2[j] =      ene[j] - nu*((f03[j+1]-f03[j]) + 0.5*tau*(D10*(g03[j+1]-g03[j])+D11*(g13[j+1]-g13[j])));

      u2[j] = mom2[j] / rho2[j];
      p2[j] = (ene2[j] - 0.5*mom2[j]*u2[j])*(gamma-1.0);

      if(p[vk1][j] < 0.0)
        printf("    (%d,%d), %g\n", k, j, p[vk1][j]);
    }

    running_info[1] = T;
    running_info[2] = 2.0;  // not half
    //WENO_5(running_info, m, h, eps, alp2, gamma, rho2, mom2, ene2, rhoI2, uI2, pI2, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R, runhist->current->trouble0);
    WENO_50(running_info, m, h, eps, alp2, gamma, rho2, mom2, ene2, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


    for(j = 0; j < m+1; ++j)
    {
      linear_GRP_solver(wave_speed, D, U, 0.0, gamma, eps,
			rho_L[j], u_L[j], 0.0, p_L[j],
			rho_R[j], u_R[j], 0.0, p_R[j],
			D_rho_L[j], D_u_L[j], 0.0, D_p_L[j],
			D_rho_R[j], D_u_R[j], 0.0, D_p_R[j],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

      g21[j] = U[0]*D[1] + U[1]*D[0];
      g22[j] = D[0]*U[1]*U[1] + 2.0*U[0]*U[1]*D[1] + D[3];
      g23[j] = (D[3]*U[1] + U[3]*D[1])*gamma/(gamma-1.0);
      g23[j] = g13[j] + 0.5*D[0]*U[1]*U[1]*U[1] + 1.5*U[0]*U[1]*U[1]*D[1];

      F21[j] = f01[j] + 0.5*tau*(D20*g01[j] + D21*g11[j] + D22*g21[j]);
      F22[j] = f02[j] + 0.5*tau*(D20*g02[j] + D21*g12[j] + D22*g22[j]);
      F23[j] = f03[j] + 0.5*tau*(D20*g03[j] + D21*g13[j] + D22*g23[j]);

      //rhoI[j] = rhoI[j] + tau*D[0];
      //momI[j] = momI[j] + tau*(D[1]*U[0] + D[0]*U[1]);
      //eneI[j] = eneI[j] + tau*(0.5*D[0]*U[1]*U[1] + U[0]*U[1]*D[1] + D[3]/(gamma-1.0));

      //uI[j] = uI[j] + tau*D[1];
      //pI[j] = pI[j] + tau*D[3];
    }
    for(j = 0; j < m; ++j)
    {
      rho[vk1][j] = rho[vk0][j] - nu*(F21[j+1]-F21[j]);
           mom[j] =      mom[j] - nu*(F22[j+1]-F22[j]);
           ene[j] =      ene[j] - nu*(F23[j+1]-F23[j]);

      u[vk1][j] = mom[j] / rho[vk1][j];
      p[vk1][j] = (ene[j] - 0.5*mom[j]*u[vk1][j])*(gamma-1.0);

      if(p[vk1][j] < 0.0)
        printf("    (%d,%d), %g\n", k, j, p[vk1][j]);
    }


    running_info[1] = T;
    running_info[2] = 0.0;  // not half
    //WENO_5(running_info, m, h, eps, alp2, gamma, rho[vk1], mom, ene, rhoI, uI, pI, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R, runhist->current->trouble0);
    WENO_50(running_info, m, h, eps, alp2, gamma, rho[vk1], mom, ene, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);

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