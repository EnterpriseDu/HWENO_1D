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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "Riemann_solver.h"

#ifndef L_STR
#include "file_io_local.h"
#endif

#include "reconstruction.h"



int LF4_WENO5_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label)
{
  delete_runHist(runhist);
  double const gamma     = OPT.gamma;  // the constant of the perfect gas
  double const gamma1 = gamma-1.0;
  double const CFL       = OPT.CFL;  // CFL number
  double const eps       = OPT.eps;  // the largest value could be treat as zero
  double const tol       = OPT.tol;

  double const alp1      = OPT.alp1;
  double const alp2      = OPT.alp2;
  double const bet1      = OPT.bet1;
  double const bet2      = OPT.bet2;
  double const modifier  = modifier;
  double const threshold = OPT.threshold;
  double const thickness = OPT.thickness;

  int const    MaxStp    = OPT.MaxStp;  //the number of time steps
  int const    bod       = OPT.bod;
  int const    Primative = OPT.Primative;
  int const    Deri      = OPT.Deri;
  int const    Limiter   = OPT.Limiter;
  int const    Decomp    = OPT.Decomp;
  double const *TIME     = OPT.TIME;
  int const    nTIME     = OPT.nTIME;
  int itTIME = 0, output = 0;


  int i = 0, j = 0, k = 1, it = 0;
  int state, len = 0;
  char scheme[L_STR] = "LF4W5\0";
  char version[L_STR], err_msg[L_STR];
  int SWITCHES[3];
  SWITCHES[0] = Deri;
  SWITCHES[1] = Decomp;
  SWITCHES[2] = Limiter;
  strcpy(version, add_mkdir);
  strcpy(add_mkdir, "../SOLUTION/\0");
  state = make_directory(add_mkdir, label, scheme, version, m, 1, SWITCHES, err_msg);
  if(state)
  {
    printf("%s", err_msg);
    exit(state);
  }
  FILE * fp_out;
  char add_out[L_STR+L_STR] = "\0", strIDX[L_STR];
  int it_out;
  double * DATA[OPT.nInitValue];
  int output_flag[OPT.nInitValue];
  for(j = 0; j < OPT.nInitValue; ++j)
    output_flag[j] = OPT.output_flag[j];
  output_flag[4] = 0;
  DATA[0] = rho;
  DATA[1] = u;
  DATA[2] = p;


  printf("The output time points are:\n  ");
  for(j = 0; j < nTIME; ++j)
    printf("| %g ", TIME[j]);
  printf("|\n\n");
  printf("===========================\n");
  printf("The scheme [%s] started.\n", scheme);


  clock_t tic, toc;
  int num_print_1 = 0, num_print_2 = 0, it_print;
  double current_cpu_time = 0.0, sum_cpu_time = 0.0, T = 0.0;
  
  
  double running_info[N_RUNNING];
  running_info[3] = bod;        //the boundary condition
  running_info[4] = Deri;       //the choice of the direvative reconstruction
  running_info[5] = Decomp;     //use the charactoristic decomposition or not
  running_info[6] = Limiter;    //use the limiter or not
  running_info[7] = threshold;  //threshold
  int recon_info[3];
  recon_info[1] = Decomp;  //characteristic decomposition or not


  double mom[m], ene[m];

  double rho_1[m], mom_1[m], ene_1[m];
  double rho_2[m], mom_2[m], ene_2[m];
  double rho_3[m], mom_3[m], ene_3[m];
  double rho_L[m+1], mom_L[m+1], ene_L[m+1], u_L, p_L, c_L;
  double rho_R[m+1], mom_R[m+1], ene_R[m+1], u_R, p_R, c_R;
  double VISC, visc;

  double f01[m+1], f02[m+1], f03[m+1];
  double f11[m+1], f12[m+1], f13[m+1];
  double f21[m+1], f22[m+1], f23[m+1];
  double f31[m+1], f32[m+1], f33[m+1];
  double df01[m+1], df02[m+1], df03[m+1]; // DUAL flux
  double df11[m+1], df12[m+1], df13[m+1]; // DUAL flux
  double fnv1, fnv2, fnv3, v1, v2, v3;
  // flux without viscosity, and viscosity

  double c, sigma, speed_max;  /* speed_max denote the largest character
					 * speed at each time step
					 */
  double tau, half_tau, alp, bet;


double a10 = 1.0;
double a20 = 649.0/1600.0, a21 = 951.0/1600.0;
double a30 = 53989.0/2500000.0, a31 = 4806213.0/20000000.0, a32 = 23619.0/32000.0;
double a40 = 0.2, a41 = 6127.0/30000.0, a42 = 7873.0/30000.0, a43 = 1.0/3.0;
double b10 = 0.5;
double b20 = -10890423.0/25193600.0, b21 = 5000.0/7873.0;
double b30 = -102261.0/5000000.0, b31 = -5121.0/20000.0, b32 = 7873.0/10000.0;
double b40 = 0.1, b41 = 1.0/6.0, b42 = 0.0, b43 = 1.0/6.0;



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


//------------THE MAIN LOOP-------------
  for(k = 1; k <= MaxStp; ++k)
  {
    if(T+eps > TIME[itTIME])
      if(itTIME == nTIME-1)
	break;
      else
	output = 1;


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
	speed_max = ((speed_max < sigma) ? sigma : speed_max);
      }
    tau = (CFL * h) / speed_max;
    if(T+tau > TIME[itTIME]){tau = TIME[itTIME]-T; T = TIME[itTIME];} else{T += tau;}
    half_tau = 0.5*tau;
    runhist->current->time[0] = tau;


    for(it_print = 0; it_print < num_print_1; ++it_print)
      printf(" ");
    printf("\r");
    fflush(stdout);
    num_print_1 = printf("%d | %g : %g | %g : %g", k, tau, T, current_cpu_time, sum_cpu_time);
    fflush(stdout);
    printf("\r");
    tic = clock();


    //======FIRST=====================
    running_info[1] = T - tau;
    running_info[2] = 0.0;
    WENO_5_LF(running_info, m, h, eps, alp2, gamma, rho, mom, ene, rho_L, rho_R, mom_L, mom_R, ene_L, ene_R);
    for(j = 0; j < m+1; ++j)
    {
      u_L = mom_L[j] / rho_L[j];
      p_L = (ene_L[j] - 0.5*mom_L[j]*u_L)*(gamma1);
      c_L = sqrt(gamma*p_L/rho_L[j]);
      u_R = mom_R[j] / rho_R[j];
      p_R = (ene_R[j] - 0.5*mom_R[j]*u_R)*(gamma1);
      c_R = sqrt(gamma*p_R/rho_R[j]);

      VISC = fabs(u_L + c_L);
      visc = fabs(u_L - c_L);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R - c_R);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R + c_R);
      if(visc > VISC)
	VISC = visc;

      fnv1 = mom_R[j]           + mom_L[j];
      fnv2 = mom_R[j]*u_R+p_R   + mom_L[j]*u_L+p_L;
      fnv3 = (ene_R[j]+p_R)*u_R + (ene_L[j]+p_L)*u_L;
      v1 = VISC*(rho_R[j]-rho_L[j]);
      v2 = VISC*(mom_R[j]-mom_L[j]);
      v3 = VISC*(ene_R[j]-ene_L[j]);

      f01[j] = 0.5*(fnv1 - v1); df01[j] = 0.5*(fnv1 + v1);
      f02[j] = 0.5*(fnv2 - v2); df02[j] = 0.5*(fnv2 + v2);
      f03[j] = 0.5*(fnv3 - v3); df03[j] = 0.5*(fnv3 + v3);
    }

    for(j = 0; j < m; ++j)
    {
      rho_1[j] = rho[j] - half_tau*(f01[j+1]-f01[j])/h;
      mom_1[j] = mom[j] - half_tau*(f02[j+1]-f02[j])/h;
      ene_1[j] = ene[j] - half_tau*(f03[j+1]-f03[j])/h;
    }

    //======SECOND=====================
    running_info[1] = T - 0.75*tau;
    running_info[2] = 1.0;
    WENO_5_LF(running_info, m, h, eps, alp2, gamma, rho_1, mom_1, ene_1, rho_L, rho_R, mom_L, mom_R, ene_L, ene_R);
    for(j = 0; j < m+1; ++j)
    {
      u_L = mom_L[j] / rho_L[j];
      p_L = (ene_L[j] - 0.5*mom_L[j]*u_L)*(gamma1);
      c_L = sqrt(gamma*p_L/rho_L[j]);
      u_R = mom_R[j] / rho_R[j];
      p_R = (ene_R[j] - 0.5*mom_R[j]*u_R)*(gamma1);
      c_R = sqrt(gamma*p_R/rho_R[j]);

      VISC = fabs(u_L + c_L);
      visc = fabs(u_L - c_L);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R - c_R);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R + c_R);
      if(visc > VISC)
	VISC = visc;

      fnv1 = mom_R[j]           + mom_L[j];
      fnv2 = mom_R[j]*u_R+p_R   + mom_L[j]*u_L+p_L;
      fnv3 = (ene_R[j]+p_R)*u_R + (ene_L[j]+p_L)*u_L;
      v1 = VISC*(rho_R[j]-rho_L[j]);
      v2 = VISC*(mom_R[j]-mom_L[j]);
      v3 = VISC*(ene_R[j]-ene_L[j]);

      f11[j] = 0.5*(fnv1 - v1); df11[j] = 0.5*(fnv1 + v1);
      f12[j] = 0.5*(fnv2 - v2); df12[j] = 0.5*(fnv2 + v2);
      f13[j] = 0.5*(fnv3 - v3); df13[j] = 0.5*(fnv3 + v3);
    }
    for(j = 0; j < m; ++j)
    {
      // b20 < 0
      rho_2[j] = (a20*rho[j]+a21*rho_1[j]) - tau*(b20*(df01[j+1]-df01[j])+b21*(f11[j+1]-f11[j]))/h;
      mom_2[j] = (a20*mom[j]+a21*mom_1[j]) - tau*(b20*(df02[j+1]-df02[j])+b21*(f12[j+1]-f12[j]))/h;
      ene_2[j] = (a20*ene[j]+a21*ene_1[j]) - tau*(b20*(df03[j+1]-df03[j])+b21*(f13[j+1]-f13[j]))/h;
    }

    //======THIRD=====================
    running_info[1] = T - 0.5*tau;
    running_info[2] = 2.0;
    WENO_5_LF(running_info, m, h, eps, alp2, gamma, rho_2, mom_2, ene_2, rho_L, rho_R, mom_L, mom_R, ene_L, ene_R);
    for(j = 0; j < m+1; ++j)
    {
      u_L = mom_L[j] / rho_L[j];
      p_L = (ene_L[j] - 0.5*mom_L[j]*u_L)*(gamma1);
      c_L = sqrt(gamma*p_L/rho_L[j]);
      u_R = mom_R[j] / rho_R[j];
      p_R = (ene_R[j] - 0.5*mom_R[j]*u_R)*(gamma1);
      c_R = sqrt(gamma*p_R/rho_R[j]);

      VISC = fabs(u_L + c_L);
      visc = fabs(u_L - c_L);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R - c_R);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R + c_R);
      if(visc > VISC)
	VISC = visc;

      f21[j] = 0.5*(mom_R[j]           + mom_L[j]           - VISC*(rho_R[j]-rho_L[j]));
      f22[j] = 0.5*(mom_R[j]*u_R+p_R   + mom_L[j]*u_L+p_L   - VISC*(mom_R[j]-mom_L[j]));
      f23[j] = 0.5*((ene_R[j]+p_R)*u_R + (ene_L[j]+p_L)*u_L - VISC*(ene_R[j]-ene_L[j]));
    }
    for(j = 0; j < m; ++j)
    {
      // b30 < 0, b31 < 0
      rho_3[j] = (a30*rho[j]+a31*rho_1[j]+a32*rho_2[j]) - tau*(b30*(df01[j+1]-df01[j])+b31*(df11[j+1]-df11[j])+b32*(f21[j+1]-f21[j]))/h;
      mom_3[j] = (a30*mom[j]+a31*mom_1[j]+a32*mom_2[j]) - tau*(b30*(df02[j+1]-df02[j])+b31*(df12[j+1]-df12[j])+b32*(f22[j+1]-f22[j]))/h;
      ene_3[j] = (a30*ene[j]+a31*ene_1[j]+a32*ene_2[j]) - tau*(b30*(df03[j+1]-df03[j])+b31*(df13[j+1]-df13[j])+b32*(f23[j+1]-f23[j]))/h;
    }

    //======FORTH=====================
    running_info[1] = T - 0.25*tau;
    running_info[2] = 3.0;
    WENO_5_LF(running_info, m, h, eps, alp2, gamma, rho_3, mom_3, ene_3, rho_L, rho_R, mom_L, mom_R, ene_L, ene_R);
    for(j = 0; j < m+1; ++j)
    {
      u_L = mom_L[j] / rho_L[j];
      p_L = (ene_L[j] - 0.5*mom_L[j]*u_L)*(gamma1);
      c_L = sqrt(gamma*p_L/rho_L[j]);
      u_R = mom_R[j] / rho_R[j];
      p_R = (ene_R[j] - 0.5*mom_R[j]*u_R)*(gamma1);
      c_R = sqrt(gamma*p_R/rho_R[j]);

      VISC = fabs(u_L + c_L);
      visc = fabs(u_L - c_L);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R - c_R);
      if(visc > VISC)
	VISC = visc;
      visc = fabs(u_R + c_R);
      if(visc > VISC)
	VISC = visc;

      f31[j] = 0.5*(mom_R[j]           + mom_L[j]           - VISC*(rho_R[j]-rho_L[j]));
      f32[j] = 0.5*(mom_R[j]*u_R+p_R   + mom_L[j]*u_L+p_L   - VISC*(mom_R[j]-mom_L[j]));
      f33[j] = 0.5*((ene_R[j]+p_R)*u_R + (ene_L[j]+p_L)*u_L - VISC*(ene_R[j]-ene_L[j]));
    }


//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
      rho[j] = (a40*rho[j]+a41*rho_1[j]+a42*rho_2[j]+a43*rho_3[j]) - tau*(b40*(f01[j+1]-f01[j])+b41*(f11[j+1]-f11[j])+b42*(f21[j+1]-f21[j])+b43*(f31[j+1]-f31[j]))/h;
      mom[j] = (a40*mom[j]+a41*mom_1[j]+a42*mom_2[j]+a43*mom_3[j]) - tau*(b40*(f02[j+1]-f02[j])+b41*(f12[j+1]-f12[j])+b42*(f22[j+1]-f22[j])+b43*(f32[j+1]-f32[j]))/h;
      ene[j] = (a40*ene[j]+a41*ene_1[j]+a42*ene_2[j]+a43*ene_3[j]) - tau*(b40*(f03[j+1]-f03[j])+b41*(f13[j+1]-f13[j])+b42*(f23[j+1]-f23[j])+b43*(f33[j+1]-f33[j]))/h;

      u[j] = mom[j] / rho[j];
      p[j] = (ene[j] - 0.5*mom[j]*u[j])*(gamma-1.0);
    }

    toc = clock();
    if(output)
    {
      printf("\nOutput at time (%d,%g).\n", itTIME, TIME[itTIME]);
      for(it_out = 0; it_out < OPT.nInitValue; ++it_out)
      {
	if(!(output_flag[it_out]))
	  continue;

	strcpy(add_out, add_mkdir);
	strcat(add_out, OPT.addInitValue[it_out]);
	sprintf(strIDX, "_%04d", itTIME);
	strcat(add_out, strIDX);
	strcat(add_out, ".txt\0");

	if((fp_out = fopen(add_out, "w")) == 0)
	{
	  sprintf(err_msg, "Cannot open solution output file: %s!\n", add_out);
	  exit(999);
	}
	for(j = 0; j < OPT.sizeInitValue[it_out]; ++j)
	  fprintf(fp_out, "%.18lf\t", DATA[it_out][j]);

	fclose(fp_out);
      }

      output = 0;
      ++itTIME;
    }
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
