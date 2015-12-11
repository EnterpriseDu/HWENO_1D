#include <math.h>
#include <stdio.h>
#include "reconstruction.h"




void HWENO_5
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double const rhoI[], double const momI[], double const eneI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  int const decmop = (int)running_info[6];
  double const threshold = running_info[7];

  double Q1[6], Q2[6], Q3[6], QI1[5], QI2[5], QI3[5], DQ1[6], DQ2[6], DQ3[6];
  double SL, SR, S0, pn1, p0, pp1;

  double H_star, u_star, c_star, u, p, H[m+4];

  double W1[m+4], W2[m+4], W3[m+4];
  double WI1[m+5], WI2[m+5], WI3[m+5];
  for(j = 0; j < m; ++j)
  {
    W1[j+2] = rho[j]; W2[j+2] = mom[j]; W3[j+2] = ene[j];
    WI1[j+2] = rhoI[j]; WI2[j+2] = momI[j]; WI3[j+2] = eneI[j];
  }
  WI1[m+2] = rhoI[m]; WI2[m+2] = momI[m]; WI3[m+2] = eneI[m];

  if(bod < 0)
  {
    W1[0]   = rho[1];  W1[1]   = rho[0];
    W2[0]   =-mom[1];  W2[1]   =-mom[0];
    W3[0]   = ene[1];  W3[1]   = ene[0];
    W1[m+2] = rho[m-1];W1[m+3] = rho[m-2];
    W2[m+2] =-mom[m-1];W2[m+3] =-mom[m-2];
    W3[m+2] = ene[m-1];W3[m+3] = ene[m-2];
    WI1[0]   = rhoI[2];  WI1[1]   = rhoI[1];
    WI2[0]   =-momI[2];  WI2[1]   =-momI[1];
    WI3[0]   = eneI[2];  WI3[1]   = eneI[1];
    WI1[m+3] = rhoI[m-1];WI1[m+4] = rhoI[m-2];
    WI2[m+3] =-momI[m-1];WI2[m+4] =-momI[m-2];
    WI3[m+3] = eneI[m-1];WI3[m+4] = eneI[m-2];
    WI2[2] = 0.0; WI2[m+2]=0.0;
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];
    W1[m+2] = rho[m-1];W1[m+3] = rho[m-1];
    W2[m+2] = mom[m-1];W2[m+3] = mom[m-1];
    W3[m+2] = ene[m-1];W3[m+3] = ene[m-1];
    WI1[0]   = rhoI[0];  WI1[1]   = rhoI[0];
    WI2[0]   = momI[0];  WI2[1]   = momI[0];
    WI3[0]   = eneI[0];  WI3[1]   = eneI[0];
    WI1[m+3] = rhoI[m];  WI1[m+4] = rhoI[m];
    WI2[m+3] = momI[m];  WI2[m+4] = momI[m];
    WI3[m+3] = eneI[m];  WI3[m+4] = eneI[m];
  }
  else
  {
    W1[0]   = rho[m-2]; W1[1]   = rho[m-1];
    W2[0]   = mom[m-2]; W2[1]   = mom[m-1];
    W3[0]   = ene[m-2]; W3[1]   = ene[m-1];
    W1[m+2] = rho[0];   W1[m+3] = rho[1];
    W2[m+2] = mom[0];   W2[m+3] = mom[1];
    W3[m+2] = ene[0];   W3[m+3] = ene[1];
    WI1[0]  = rhoI[m-2];WI1[1]  = rhoI[m-1];
    WI2[0]  = momI[m-2];WI2[1]  = momI[m-1];
    WI3[0]  = eneI[m-2];WI3[1]  = eneI[m-1];
    WI1[m+3] = rhoI[1]; WI1[m+4] = rhoI[2];
    WI2[m+3] = momI[1]; WI2[m+4] = momI[2];
    WI3[m+3] = eneI[1]; WI3[m+4] = eneI[2];
  }

  for(j = 0; j < m+4; ++j)
  {
    u = W2[j]/W1[j];
    p = (W3[j] - 0.5*W2[j]*u)*(gamma-1.0);
    H[j] = 0.5*u*u + gamma*p/W1[j]/(gamma-1.0);
  }

  for(j = 2; j < m+3; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 5; ++k)
    {
      QI1[k] = WI1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      QI1[k]+= WI2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      QI1[k]+= WI3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      QI1[k] = QI1[k] / c_star;
      QI3[k] = WI1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      QI3[k]+= WI2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      QI3[k]+= WI3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      QI3[k] = QI3[k] / c_star;
      QI2[k] = (gamma-1.0) * ((WI2[j-2+k]-0.5*WI1[j-2+k]*u_star)*u_star - WI3[j-2+k]) / c_star/c_star + WI1[j-2+k];
    }
    for(k = 0; k < 4; ++k)
    {
      Q1[k] = W1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-2+k]-0.5*W1[j-2+k]*u_star)*u_star - W3[j-2+k]) / c_star/c_star + W1[j-2+k];

      DQ1[k] = QI1[k+1] - QI1[k];
      DQ2[k] = QI2[k+1] - QI2[k];
      DQ3[k] = QI3[k+1] - QI3[k];
    }

    if(WENOD)
    {
      local_HWENO_5_inter_d(h, Q1, DQ1);
      local_HWENO_5_inter_d(h, Q2, DQ2);
      local_HWENO_5_inter_d(h, Q3, DQ3);
    }
    else
    {
      local_HWENO_5_inter(h, Q1, DQ1);
      local_HWENO_5_inter(h, Q2, DQ2);
      local_HWENO_5_inter(h, Q3, DQ3);
      DQ1[4] = (15*(Q1[2] - Q1[1]) + Q1[0] - Q1[3]) / (12*h);
      DQ2[4] = (15*(Q2[2] - Q2[1]) + Q2[0] - Q2[3]) / (12*h);
      DQ3[4] = (15*(Q3[2] - Q3[1]) + Q3[0] - Q3[3]) / (12*h);
      DQ1[5] = DQ1[4];
      DQ2[5] = DQ2[4];
      DQ3[5] = DQ3[4];
    }

  //=====Recomposition========
    rho_L[j-2] = Q1[4] + Q2[4] + Q3[4];
    rho_R[j-2] = Q1[5] + Q2[5] + Q3[5];
    D_rho_L[j-2] = DQ1[4] + DQ2[4] + DQ3[4];
    D_rho_R[j-2] = DQ1[5] + DQ2[5] + DQ3[5];

    u_L[j-2] = u_star*rho_L[j-2] + c_star*(Q3[4] - Q1[4]);
    u_R[j-2] = u_star*rho_R[j-2] + c_star*(Q3[5] - Q1[5]);
    D_u_L[j-2] = u_star*D_rho_L[j-2] + c_star*(DQ3[4] - DQ1[4]);
    D_u_R[j-2] = u_star*D_rho_R[j-2] + c_star*(DQ3[5] - DQ1[5]);

    p_L[j-2] = H_star*(Q1[4]+Q3[4]) + u_star*c_star*(Q3[4]-Q1[4]) + 0.5*u_star*u_star*Q2[4];
    p_R[j-2] = H_star*(Q1[5]+Q3[5]) + u_star*c_star*(Q3[5]-Q1[5]) + 0.5*u_star*u_star*Q2[5];
    D_p_L[j-2] = H_star*(DQ1[4]+DQ3[4]) + u_star*c_star*(DQ3[4]-DQ1[4]) + 0.5*u_star*u_star*DQ2[4];
    D_p_R[j-2] = H_star*(DQ1[5]+DQ3[5]) + u_star*c_star*(DQ3[5]-DQ1[5]) + 0.5*u_star*u_star*DQ2[5];


    u_L[j-2] = u_L[j-2] / rho_L[j-2];
    u_R[j-2] = u_R[j-2] / rho_R[j-2];
    p_L[j-2] = (p_L[j-2] - 0.5*u_L[j-2]*u_L[j-2]*rho_L[j-2])*(gamma-1.0);
    p_R[j-2] = (p_R[j-2] - 0.5*u_R[j-2]*u_R[j-2]*rho_R[j-2])*(gamma-1.0);

    D_u_L[j-2] = (D_u_L[j-2] - u_L[j-2]*D_rho_L[j-2]) / rho_L[j-2];
    D_u_R[j-2] = (D_u_R[j-2] - u_R[j-2]*D_rho_R[j-2]) / rho_R[j-2];
    D_p_L[j-2] = (D_p_L[j-2] - 0.5*D_rho_L[j-2]*u_L[j-2]*u_L[j-2] - rho_L[j-2]*u_L[j-2]*D_u_L[j-2])*(gamma-1.0);
    D_p_R[j-2] = (D_p_R[j-2] - 0.5*D_rho_R[j-2]*u_R[j-2]*u_R[j-2] - rho_R[j-2]*u_R[j-2]*D_u_R[j-2])*(gamma-1.0);
  }
}




void HWENO_5_limited
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double const rhoI[], double const momI[], double const eneI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  double const threshold = running_info[6];

  double Q1[6], Q2[6], Q3[6], QI1[5], QI2[5], QI3[5], DQ1[6], DQ2[6], DQ3[6];
  double u[m+4], p[m+4], H[m+4], Entp[m+4], EntpI[m+5];
  double EntpL, EntpR;
  double deltaP, DP, theta = 25.0;
  double Qn1, Q0, Qp1, SL, SR, S0;
  int flag;

  double H_star, u_star, c_star;

  double W1[m+4], W2[m+4], W3[m+4];
  double WI1[m+5], WI2[m+5], WI3[m+5];
  for(j = 0; j < m; ++j)
  {
    W1[j+2] = rho[j]; W2[j+2] = mom[j]; W3[j+2] = ene[j];
    trouble[j] = 0;
  }
  for(j = 0; j < m+1; ++j)
  {
    WI1[j+2] = rhoI[j]; WI2[j+2] = momI[j]; WI3[j+2] = eneI[j];
    EntpI[j+2] = pI[j] / pow(rhoI[j], gamma);
  }

  if(bod < 0)
  {
    W1[0]   = rho[1];  W1[1]   = rho[0];
    W2[0]   =-mom[1];  W2[1]   =-mom[0];
    W3[0]   = ene[1];  W3[1]   = ene[0];
    W1[m+2] = rho[m-1];W1[m+3] = rho[m-2];
    W2[m+2] =-mom[m-1];W2[m+3] =-mom[m-2];
    W3[m+2] = ene[m-1];W3[m+3] = ene[m-2];
    WI1[0]   = rhoI[2];  WI1[1]   = rhoI[1];
    WI2[0]   =-momI[2];  WI2[1]   =-momI[1];
    WI3[0]   = eneI[2];  WI3[1]   = eneI[1];
    WI1[m+3] = rhoI[m-1];WI1[m+4] = rhoI[m-2];
    WI2[m+3] =-momI[m-1];WI2[m+4] =-momI[m-2];
    WI3[m+3] = eneI[m-1];WI3[m+4] = eneI[m-2];
    WI2[2] = 0.0; WI2[m+2]=0.0;
    EntpI[0]   = EntpI[4];  EntpI[1]   = EntpI[3];
    EntpI[m+3]   = EntpI[m+1];  EntpI[m+4]   = EntpI[m];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];
    W1[m+2] = rho[m-1];W1[m+3] = rho[m-1];
    W2[m+2] = mom[m-1];W2[m+3] = mom[m-1];
    W3[m+2] = ene[m-1];W3[m+3] = ene[m-1];
    WI1[0]   = rhoI[0];  WI1[1]   = rhoI[0];
    WI2[0]   = momI[0];  WI2[1]   = momI[0];
    WI3[0]   = eneI[0];  WI3[1]   = eneI[0];
    WI1[m+3] = rhoI[m];  WI1[m+4] = rhoI[m];
    WI2[m+3] = momI[m];  WI2[m+4] = momI[m];
    WI3[m+3] = eneI[m];  WI3[m+4] = eneI[m];
    EntpI[0]   = EntpI[2];  EntpI[1]   = EntpI[2];
    EntpI[m+3]   = EntpI[m+2];  EntpI[m+4]   = EntpI[m+2];
  }
  else
  {
    W1[0]   = rho[m-2]; W1[1]   = rho[m-1];
    W2[0]   = mom[m-2]; W2[1]   = mom[m-1];
    W3[0]   = ene[m-2]; W3[1]   = ene[m-1];
    W1[m+2] = rho[0];   W1[m+3] = rho[1];
    W2[m+2] = mom[0];   W2[m+3] = mom[1];
    W3[m+2] = ene[0];   W3[m+3] = ene[1];
    WI1[0]  = rhoI[m-2];WI1[1]  = rhoI[m-1];
    WI2[0]  = momI[m-2];WI2[1]  = momI[m-1];
    WI3[0]  = eneI[m-2];WI3[1]  = eneI[m-1];
    WI1[m+3] = rhoI[1]; WI1[m+4] = rhoI[2];
    WI2[m+3] = momI[1]; WI2[m+4] = momI[2];
    WI3[m+3] = eneI[1]; WI3[m+4] = eneI[2];
    EntpI[0]   = EntpI[m];  EntpI[1]   = EntpI[m+1];
    EntpI[m+3]   = EntpI[3];  EntpI[m+4]   = EntpI[4];
  }

  /*
  W1[m+2] = 0.2*(cos(25.0+5.0*h) - cos(25.0));
  W1[m+2] = 1.0 + 0.2*W1[m+2];
  W1[m+3] = 0.2*(cos(25.0+10.0*h) - cos(25.0+5.0*h));
  W1[m+3] = 1.0 + 0.2*W1[m+3];
  WI1[m+2] = 1.0 + 0.2*sin(25.0);
  WI1[m+3] = 1.0 + 0.2*sin(25.0+5.0*h);
  WI1[m+4] = 1.0 + 0.2*sin(25.0+10.0*h);
  EntpI[m+2] = 1.0 / pow(WI1[m+2], gamma);
  EntpI[m+3] = 1.0 / pow(WI1[m+3], gamma);
  EntpI[m+4] = 1.0 / pow(WI1[m+4], gamma);
  */


  for(j = 0; j < m+4; ++j)
  {
    u[j] = W2[j]/W1[j];
    p[j] = (W3[j] - 0.5*W2[j]*u[j])*(gamma-1.0);
    H[j] = 0.5*u[j]*u[j] + gamma*p[j]/W1[j]/(gamma-1.0);
    Entp[j] = p[j] / pow(W1[j], gamma);
  }

  for(j = 2; j < m+3; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 5; ++k)
    {
      QI1[k] = WI1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      QI1[k]+= WI2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      QI1[k]+= WI3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      QI1[k] = QI1[k] / c_star;
      QI3[k] = WI1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      QI3[k]+= WI2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      QI3[k]+= WI3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      QI3[k] = QI3[k] / c_star;
      QI2[k] = (gamma-1.0) * ((WI2[j-2+k]-0.5*WI1[j-2+k]*u_star)*u_star - WI3[j-2+k]) / c_star/c_star + WI1[j-2+k];
    }
    for(k = 0; k < 4; ++k)
    {
      Q1[k] = W1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-2+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-2+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-2+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-2+k]-0.5*W1[j-2+k]*u_star)*u_star - W3[j-2+k]) / c_star/c_star + W1[j-2+k];

      DQ1[k] = QI1[k+1] - QI1[k];
      DQ2[k] = QI2[k+1] - QI2[k];
      DQ3[k] = QI3[k+1] - QI3[k];
    }
    if(WENOD)
    {
      local_HWENO_5_inter_d(h, Q1, DQ1);
      local_HWENO_5_inter_d(h, Q2, DQ2);
      local_HWENO_5_inter_d(h, Q3, DQ3);
    }
    else
    {
      local_HWENO_5_inter(h, Q1, DQ1);
      local_HWENO_5_inter(h, Q2, DQ2);
      local_HWENO_5_inter(h, Q3, DQ3);
      DQ1[4] = (15*(Q1[2] - Q1[1]) + Q1[0] - Q1[3]) / (12*h);
      DQ2[4] = (15*(Q2[2] - Q2[1]) + Q2[0] - Q2[3]) / (12*h);
      DQ3[4] = (15*(Q3[2] - Q3[1]) + Q3[0] - Q3[3]) / (12*h);
      DQ1[5] = DQ1[4];
      DQ2[5] = DQ2[4];
      DQ3[5] = DQ3[4];
    }

  //=====Recomposition========
    rho_L[j-2] = Q1[4] + Q2[4] + Q3[4];
    rho_R[j-2] = Q1[5] + Q2[5] + Q3[5];
    D_rho_L[j-2] = DQ1[4] + DQ2[4] + DQ3[4];
    D_rho_R[j-2] = DQ1[5] + DQ2[5] + DQ3[5];

    u_L[j-2] = u_star*rho_L[j-2] + c_star*(Q3[4] - Q1[4]);
    u_R[j-2] = u_star*rho_R[j-2] + c_star*(Q3[5] - Q1[5]);
    D_u_L[j-2] = u_star*D_rho_L[j-2] + c_star*(DQ3[4] - DQ1[4]);
    D_u_R[j-2] = u_star*D_rho_R[j-2] + c_star*(DQ3[5] - DQ1[5]);

    p_L[j-2] = H_star*(Q1[4]+Q3[4]) + u_star*c_star*(Q3[4]-Q1[4]) + 0.5*u_star*u_star*Q2[4];
    p_R[j-2] = H_star*(Q1[5]+Q3[5]) + u_star*c_star*(Q3[5]-Q1[5]) + 0.5*u_star*u_star*Q2[5];
    D_p_L[j-2] = H_star*(DQ1[4]+DQ3[4]) + u_star*c_star*(DQ3[4]-DQ1[4]) + 0.5*u_star*u_star*DQ2[4];
    D_p_R[j-2] = H_star*(DQ1[5]+DQ3[5]) + u_star*c_star*(DQ3[5]-DQ1[5]) + 0.5*u_star*u_star*DQ2[5];


    u_L[j-2] = u_L[j-2] / rho_L[j-2];
    u_R[j-2] = u_R[j-2] / rho_R[j-2];
    p_L[j-2] = (p_L[j-2] - 0.5*u_L[j-2]*u_L[j-2]*rho_L[j-2])*(gamma-1.0);
    p_R[j-2] = (p_R[j-2] - 0.5*u_R[j-2]*u_R[j-2]*rho_R[j-2])*(gamma-1.0);

    D_u_L[j-2] = (D_u_L[j-2] - u_L[j-2]*D_rho_L[j-2]) / rho_L[j-2];
    D_u_R[j-2] = (D_u_R[j-2] - u_R[j-2]*D_rho_R[j-2]) / rho_R[j-2];
    D_p_L[j-2] = (D_p_L[j-2] - 0.5*D_rho_L[j-2]*u_L[j-2]*u_L[j-2] - rho_L[j-2]*u_L[j-2]*D_u_L[j-2])*(gamma-1.0);
    D_p_R[j-2] = (D_p_R[j-2] - 0.5*D_rho_R[j-2]*u_R[j-2]*u_R[j-2] - rho_R[j-2]*u_R[j-2]*D_u_R[j-2])*(gamma-1.0);



    flag = 0;
    /*
    deltaP = threshold*fabs(Q1[4] - Q1[5]);
    DP = 1e-5+fabs(30.0*(Q1[2]-Q1[1]) + (520.0/27.0)*(Q1[0]-Q1[3]) + (20.0/3.0)*(QI1[4]-QI1[3]+QI1[1]-QI1[0]));
    //if(deltaP > DP) flag = 1*limiter;

    deltaP = threshold*fabs(Q2[4] - Q2[5]);
    DP = 1e-5+fabs(30.0*(Q2[2]-Q2[1]) + (520.0/27.0)*(Q2[0]-Q2[3]) + (20.0/3.0)*(QI2[4]-QI2[3]+QI2[1]-QI2[0]));
    //if(deltaP > DP) flag = 2*limiter;

    deltaP = threshold*fabs(Q3[4] - Q3[5]);
    DP = 1e-5+fabs(30.0*(Q3[2]-Q3[1]) + (520.0/27.0)*(Q3[0]-Q3[3]) + (20.0/3.0)*(QI3[4]-QI3[3]+QI3[1]-QI3[0]));
    //if(deltaP > DP) flag = 3*limiter;
    //*/
    deltaP = threshold*fabs(rho_L[j-2] - rho_R[j-2]);
    DP = 1e-5+fabs(30.0*(W1[j]-W1[j-1]) + (520.0/27.0)*(W1[j-2]-W1[j+1]) + (20.0/3.0)*(WI1[j+2]-WI1[j+1]+WI1[j-1]-WI1[j-2]));
    if(deltaP > DP) flag = 4*limiter;
    /*
    for(k = 0; k < 5; ++k)
      QI1[k] = EntpI[j-2+k];
    for(k = 0; k < 4; ++k)
    {
      Q1[k] = Entp[j-2+k];
      DQ1[k] = QI1[k+1] - QI1[k];
    }
    local_HWENO_5_inter(h, Q1, DQ1);
    deltaP = threshold*fabs(Q1[4] - Q1[5]);
    DP = 1e-5+fabs(30.0*(Q1[2]-Q1[1]) + (520.0/27.0)*(Q1[0]-Q1[3]) + (20.0/3.0)*(QI1[4]-QI1[3]+QI1[1]-QI1[0]));
    if(deltaP > DP) flag = 5*limiter;
    //EntpL = p_L[j-2] / pow(rho_L[j-2], gamma);
    //EntpR = p_R[j-2] / pow(rho_R[j-2], gamma);
    //deltaP = threshold*fabs(EntpL[j-2] - EntpR[j-2]);
    //DP = 1e-5+fabs(30.0*(Entp[j]-Entp[j-1]) + (130.0/9.0)*(Entp[j-2]-Entp[j+1]) + (20.0/3.0)*(EntpI[j+2]-EntpI[j+1]+EntpI[j-1]-EntpI[j-2]));
    //if(deltaP > DP) trouble[j-2] = 6*limiter;
    //*/

    if(j-2)
      trouble[j-3] += flag;
    if(m+2-j)
      trouble[j-2] += flag;
  }

  for(j = 0; j < m; ++j)
  {
    if(trouble[j])
    {
      Qn1 = W1[j+1];
      Q0  = W1[j+2];
      Qp1 = W1[j+3];

      SL = Q0 - Qn1; SR = Qp1 - Q0; S0 = rhoI[j+1] - rhoI[j];
      if(SL*SR < 0.0)
      {
	D_rho_R[j] = 0.0;
	rho_R[j] = Q0;
	rho_L[j+1] = Q0;
	D_rho_L[j+1] = D_rho_R[j];
      }
      else if(SL*S0 < 0.0)
      {
	D_rho_R[j] = 0.0;
	rho_R[j] = Q0;
	rho_L[j+1] = Q0;
	D_rho_L[j+1] = D_rho_R[j];
      }
      else
      {
	SL = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	if(fabs(SL) < fabs(S0))
	  S0 = SL;
	rho_R[j] = Q0 - 0.5*S0;
	D_rho_R[j] = S0 / h;
	rho_L[j+1] = Q0 + 0.5*S0;
	D_rho_L[j+1] = D_rho_R[j];
      }

      Qn1 = u[j+1];
      Q0  = u[j+2];
      Qp1 = u[j+3];
      SL = Q0 - Qn1; SR = Qp1 - Q0; S0 = uI[j+1] - uI[j];
      if(SL*SR < 0.0)
      {
	D_u_R[j] = 0.0;
	u_R[j] = Q0;
	u_L[j+1] = Q0;
	D_u_L[j+1] = D_u_R[j];
      }
      else if(SL*S0 < 0.0)
      {
	D_u_R[j] = 0.0;
	u_R[j] = Q0;
	u_L[j+1] = Q0;
	D_u_L[j+1] = D_u_R[j];
      }
      else
      {
	SL = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	if(fabs(SL) < fabs(S0))
	  S0 = SL;
	u_R[j] = Q0 - 0.5*S0;
	D_u_R[j] = S0 / h;
	u_L[j+1] = Q0 + 0.5*S0;
	D_u_L[j+1] = D_u_R[j];
      }

      Qn1 = p[j+1];
      Q0  = p[j+2];
      Qp1 = p[j+3];
      SL = Q0 - Qn1; SR = Qp1 - Q0; S0 = pI[j+1] - pI[j];
      if(SL*SR < 0.0)
      {
	D_p_R[j] = 0.0;
	p_R[j] = Q0;
	p_L[j+1] = Q0;
	D_p_L[j+1] = D_p_R[j];
      }
      else if(SL*S0 < 0.0)
      {
	D_p_R[j] = 0.0;
	p_R[j] = Q0;
	p_L[j+1] = Q0;
	D_p_L[j+1] = D_p_R[j];
      }
      else
      {
	SL = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	if(fabs(SL) < fabs(S0))
	  S0 = SL;
	p_R[j] = Q0 - 0.5*S0;
	D_p_R[j] = S0 / h;
	p_L[j+1] = Q0 + 0.5*S0;
	D_p_L[j+1] = D_p_R[j];
      }
    }
  }
}


