#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "Riemann_solver.h"

#ifndef L_STR
#include "file_io_local.h"
#endif

#include "reconstruction.h"




int GRP2_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label)
{
  delete_runHist(runhist);
  double const gamma     = OPT.gamma;  // the constant of the perfect gas
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
  char scheme[L_STR] = "G2L2\0";
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

  int trouble0[m];
  FILE * fp_tr0;
  char add_tr0[L_STR+L_STR];
  strcpy(add_tr0, add_mkdir);
  strcat(add_tr0, "trouble0.txt\0");
  state = open_fruncate(err_msg, add_tr0, &fp_tr0);
  if(state)
  { printf("%s", err_msg);
    exit(state); }


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


  double c, stmp, D[5], U[5], wave_speed[2];
  RSboundary wL, wR;
  RSparameters RSpara;


  double mom[m], ene[m];
  double rho_L[m+1], rho_R[m+1], u_L[m+1], u_R[m+1], p_L[m+1], p_R[m+1];
  double D_rho_L[m+1], D_rho_R[m+1], D_u_L[m+1], D_u_R[m+1], D_p_L[m+1], D_p_R[m+1];

  double F1[m+1], F2[m+1], F3[m+1], rhoI[m+1], uI[m+1], pI[m+1];

  double sigma, speed_max;  /* speed_max denote the largest character
					 * speed at each time step
					 */
  double tau, half_tau, alp, bet, nu;


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
  GRP_minmod0(running_info, m, h, alp2, rho, u, p, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);
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
      if(speed_max < sigma)
        speed_max = sigma;
    }
    tau = (CFL * h) / speed_max;
    if(T+tau > TIME[itTIME]){tau = TIME[itTIME]-T; T = TIME[itTIME];} else{T += tau;}
    half_tau = 0.5*tau;
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
      wL.gamma = gamma;
      wL.rho = rho_L[j];
      wL.u = u_L[j];
      wL.v = 0.0;
      wL.p = p_L[j];
      wL.d_gamma = 0.0;
      wL.d_rho = D_rho_L[j];
      wL.d_u = D_u_L[j];
      wL.d_v = 0.0;
      wL.d_p = D_p_L[j];
      wL.t_gamma = 0.0;
      wL.t_rho = 0.0;
      wL.t_u = 0.0;
      wL.t_v = 0.0;
      wL.t_p = 0.0;
      wR.gamma = gamma;
      wR.rho = rho_R[j];
      wR.u = u_R[j];
      wR.v = 0.0;
      wR.p = p_R[j];
      wR.d_gamma = 0.0;
      wR.d_rho = D_rho_R[j];
      wR.d_u = D_u_R[j];
      wR.d_v = 0.0;
      wR.d_p = D_p_R[j];
      wR.t_gamma = 0.0;
      wR.t_rho = 0.0;
      wR.t_u = 0.0;
      wR.t_v = 0.0;
      wR.t_p = 0.0;
      RSpara.eps = eps;
      RSpara.tol = tol;
      RSpara.N = 500;
      RSpara.radius = 1.0;
      RSpara.nDim = 1;
      Euler_GRP_solver(wave_speed, D, U, 0.0, &wL, &wR, &RSpara);

      F1[j] = U[0]*U[1] + half_tau*(D[0]*U[1]+U[0]*D[1]);
      F2[j] = U[0]*U[1]*U[1] + half_tau*(D[0]*U[1]*U[1]+2.0*U[0]*U[1]*D[1]) + U[3]+half_tau*D[3];
      F3[j] = (U[3]*U[1]+half_tau*(D[3]*U[1]+U[3]*D[1]))*gamma/(gamma-1.0);
      F3[j] = 0.5*U[0]*U[1]*U[1]*U[1] + half_tau*(0.5*D[0]*U[1]*U[1]*U[1]+1.5*U[0]*U[1]*U[1]*D[1]) + F3[j];

      rhoI[j] = U[0] + tau*D[0];
        uI[j] = U[1] + tau*D[1];
        pI[j] = U[3] + tau*D[3];
    }

//===============THE CORE ITERATION=================
    for(j = 0; j < m; ++j)
    {
      rho[j] = rho[j] - nu*(F1[j+1]-F1[j]);
      mom[j] = mom[j] - nu*(F2[j+1]-F2[j]);
      ene[j] = ene[j] - nu*(F3[j+1]-F3[j]);

	u[j] = mom[j] / rho[j];
	p[j] = (ene[j] - 0.5*mom[j]*u[j])*(gamma-1.0);
    }

    running_info[1] = T;
    running_info[2] = 0.0;
    GRP_minmod(running_info, m, h, alp2, rho, u, p, rhoI, uI, pI, rho_L, rho_R, u_L, u_R, p_L, p_R, D_rho_L, D_rho_R, D_u_L, D_u_R, D_p_L, D_p_R);


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

    fprintf(fp_tr0, "%.10lf\t%.10lf\n", p[0], p[m-1]);
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

  fclose(fp_tr0);
  return k;
}
