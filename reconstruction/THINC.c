#include <math.h>
#include <stdio.h>


#include "reconstruction.h"


void THINC0
(double const running_info[], int const m, double const h, double const alp2,
 double const rho[], double const u[], double const p[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[])
{
  int i, j;
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];

  double u_min, u_max, u_bar, thickness, result[4];
  thickness = 1.6;

  double SL, SR, Stmp;
  double Drho[m], Du[m], Dp[m];
  double P1[m+2], P2[m+2], P3[m+2];
  for(j = 0; j < m; ++j)
  {
    P1[j+1] = rho[j];
    P2[j+1] =   u[j];
    P3[j+1] =   p[j];
  }

  printf("####%d####\n", K);

  if(bod < 0)
  {
    P1[0]   = rho[0];  P2[0]   =-u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] =-u[m-1];P3[m+1] = p[m-1];
  }
  else if(bod > 0)
  {
    P1[0]   = rho[0];  P2[0]   = u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] = u[m-1];P3[m+1] = p[m-1];
  }
  else
  {
    P1[0]   = rho[m-1];P2[0]   = u[m-1];P3[0]   = p[m-1];
    P1[m+1] = rho[0];  P2[m+1] = u[0];  P3[m+1] = p[0];
  }


  for(j = 1; j < m+1; ++j)
  {
    if(!trouble[j-1])
      continue;

    SR = (P1[j+1] - P1[j]);
    SL = (P1[j] - P1[j-1]);
    if(SR*SL < 0.0)
    {
      rho_R[j-1] = P1[j];
      rho_L[j] = P1[j];
      D_rho_R[j-1] = 0.0;
      D_rho_L[j] = 0.0;
    }
    else if(SL > 0.0)
    {
      u_min = P1[j-1];
      u_bar = P1[j];
      u_max = P1[j+1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      rho_R[j-1] = result[0];
      rho_L[j] = result[1];
      D_rho_R[j-1] = result[2];
      D_rho_L[j] = result[3];
    }
    else
    {
      u_min = P1[j+1];
      u_bar = P1[j];
      u_max = P1[j-1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      rho_R[j-1] = result[1];
      rho_L[j] = result[0];
      D_rho_R[j-1] = -result[3];
      D_rho_L[j] = -result[2];
    }

    SR = (P2[j+1] - P2[j]);
    SL = (P2[j] - P2[j-1]);
    if(SR*SL < 0.0)
    {
      u_R[j-1] = P2[j];
      u_L[j] = P2[j];
      D_u_R[j-1] = 0.0;
      D_u_L[j] = 0.0;
    }
    else if(SL > 0.0)
    {
      u_min = P2[j-1];
      u_bar = P2[j];
      u_max = P2[j+1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      u_R[j-1] = result[0];
      u_L[j] = result[1];
      D_u_R[j-1] = result[2];
      D_u_L[j] = result[3];
    }
    else
    {
      u_min = P2[j+1];
      u_bar = P2[j];
      u_max = P2[j-1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      u_R[j-1] = result[1];
      u_L[j] = result[0];
      D_u_R[j-1] = -result[3];
      D_u_L[j] = -result[2];
    }

    SR = (P3[j+1] - P3[j]);
    SL = (P3[j] - P3[j-1]);
    if(SR*SL < 0.0)
    {
      p_R[j-1] = P3[j];
      p_L[j] = P3[j];
      D_p_R[j-1] = 0.0;
      D_p_L[j] = 0.0;
    }
    else if(SL > 0.0)
    {
      u_min = P3[j-1];
      u_bar = P3[j];
      u_max = P3[j+1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      p_R[j-1] = result[0];
      p_L[j] = result[1];
      D_p_R[j-1] = result[2];
      D_p_L[j] = result[3];
    }
    else
    {
      u_min = P3[j+1];
      u_bar = P3[j];
      u_max = P3[j-1];
      THINC_local(result, u_min, u_max-u_min, u_bar, thickness, h);
      p_R[j-1] = result[1];
      p_L[j] = result[0];
      D_p_R[j-1] = -result[3];
      D_p_L[j] = -result[2];
    }
    //*
    printf("***%d***\n", j-1);
    printf("Density:\n");
    printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", P1[j-1], rho_L[j-1], rho_R[j-1], P1[j], rho_L[j], rho_R[j], P1[j+1]);
    printf("\t\t%g\t%g\t\t%g\t%g\n", D_rho_L[j-1], D_rho_R[j-1], D_rho_L[j], D_rho_R[j]);
    printf("Velocity:\n");
    printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", P2[j-1], u_L[j-1], u_R[j-1], P2[j], u_L[j], u_R[j], P2[j+1]);
    printf("\t\t%g\t%g\t\t%g\t%g\n", D_u_L[j-1], D_u_R[j-1], D_u_L[j], D_u_R[j]);
    printf("Pressure:\n");
    printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", P3[j-1], p_L[j-1], p_R[j-1], P3[j], p_L[j], p_R[j], P3[j+1]);
    printf("\t\t%g\t%g\t\t%g\t%g\n", D_p_L[j-1], D_p_R[j-1], D_p_L[j], D_p_R[j]);
    printf("\n");//*/
  }



  if(bod < 0)
  {
    if(trouble[m-1])
    {
      rho_R[m] = rho_L[m];
      u_R[m] =  -u_L[m];
      p_R[m] =   p_L[m];
      D_rho_R[m] = -D_rho_L[m];
      D_u_R[m] =    D_u_L[m];
      D_p_R[m] =   -D_p_L[m];
    }

    if(trouble[0])
    {
      rho_L[0] = rho_R[0];
      u_L[0] =  -u_R[0];
      p_L[0] =   p_R[0];
      D_rho_L[0] = -D_rho_R[0];
      D_u_L[0] =    D_u_R[0];
      D_p_L[0] =   -D_p_R[0];
    }
  }
  else if(bod > 0)
  {
    if(trouble[m-1])
    {
      rho_R[m] = rho_R[m-1];
      u_R[m] =   u_R[m-1];
      p_R[m] =   p_R[m-1];
      D_rho_R[m] = D_rho_R[m-1];
      D_u_R[m] =   D_u_R[m-1];
      D_p_R[m] =   D_p_R[m-1];
    }

    if(trouble[0])
    {
      rho_L[0] = rho_L[1];
      u_L[0] =   u_L[1];
      p_L[0] =   p_L[1];
      D_rho_L[0] = D_rho_L[1];
      D_u_L[0] =   D_u_L[1];
      D_p_L[0] =   D_p_L[1];
    }
  }
  else
  {
    if(trouble[m-1])
    {
      rho_R[m] = rho_R[0];
      u_R[m] =   u_R[0];
      p_R[m] =   p_R[0];
      D_rho_R[m] = D_rho_R[0];
      D_u_R[m] =   D_u_R[0];
      D_p_R[m] =   D_p_R[0];
    }

    if(trouble[0])
    {
      rho_L[0] = rho_L[m];
      u_L[0] =   u_L[m];
      p_L[0] =   p_L[m];
      D_rho_L[0] = D_rho_L[m];
      D_u_L[0] =   D_u_L[m];
      D_p_L[0] =   D_p_L[m];
    }
  }
}





void THINC_local(double result[], double const u_min, double const u_jump, double const u_bar, double const thickness, double const h)
{
  double eps = 1e-20;
  double B, A;
  double idx, exp_idx_0, exp_idx_1, tanh_beta, sinh_beta;


  tanh_beta = tanh(thickness);
  sinh_beta = sinh(thickness);
  idx = (u_bar-u_min+eps)/(u_jump+eps);
  //idx = thickness*(2.0*(u_bar-u_min+eps)/(u_jump+eps)-1);
  B = exp(thickness*(2.0*idx - 1.0));
  A = (B/cosh(thickness)-1.0) / tanh_beta;

  result[0] = u_min + 0.5*u_jump*(1.0+A);
  result[1] = u_min + 0.5*u_jump*(1.0+((tanh_beta+A)/(1.0+A*tanh_beta)));

  exp_idx_0 = exp(2.0*thickness*idx);
  exp_idx_1 = exp(2.0*thickness*(idx-1.0));
  result[2] = -0.5*u_jump*thickness/h/(sinh_beta*sinh_beta);
  result[3] = result[2]*(1.0/exp_idx_0 - 1.0)*(1.0/exp_idx_1 - 1.0);
  result[2] = result[2]*(exp_idx_0 - 1.0)*(exp_idx_1 - 1.0);
}
