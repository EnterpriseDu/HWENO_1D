#include <math.h>

#include "reconstruction.h"





void WENO_5_limited
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[], double const rhoI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  int const decmop = (int)running_info[6];
  double const threshold = running_info[7];


  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double PI1[m+3], PI2[m+3], PI3[m+3];
  double Q1[8], Q2[8], Q3[8], DQQ[3][2];
  double u[m+6], p[m+6], H[m+6], Entp[m+6];
  double EntpL, EntpR;
  double Qn1, Q0, Qp1, SL, SR, S0;
  double deltaP, DP;
  int flag;

  double H_star, u_star, c_star;
  for(j = 0; j < m; ++j)
  {
    W1[j+3] = rho[j];
    W2[j+3] = mom[j];
    W3[j+3] = ene[j];
    PI1[j+1] = rhoI[j]; PI2[j+1] = uI[j]; PI3[j+1] = pI[j];
    trouble[j] = 0;
  }
  PI1[m+1] = rhoI[m]; PI2[m+1] = uI[m]; PI3[m+1] = pI[m];

  if(bod < 0)
  {
    W1[0]   = rho[2];  W1[1]   = rho[1];  W1[2]   = rho[0];
    W2[0]   =-mom[2];  W2[1]   =-mom[1];  W2[2]   =-mom[0];
    W3[0]   = ene[2];  W3[1]   = ene[1];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-2];W1[m+5] = rho[m-3];
    W2[m+3] =-mom[m-1];W2[m+4] =-mom[m-2];W2[m+5] =-mom[m-3];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-2];W3[m+5] = ene[m-3];
    PI1[0] = rhoI[1]; PI1[m+2] = rhoI[m-1];
    PI2[0] =-  uI[1]; PI2[m+2] =-  uI[m-1]; PI2[1] = 0.0; PI2[m+1]=0.0;
    PI3[0] =   pI[1]; PI3[m+2] =   pI[m-1];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];  W1[2]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];  W2[2]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-1];W1[m+5] = rho[m-1];
    W2[m+3] = mom[m-1];W2[m+4] = mom[m-1];W2[m+5] = mom[m-1];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-1];W3[m+5] = ene[m-1];
    PI1[0] = rhoI[0]; PI1[m+2] = rhoI[m];
    PI2[0] =   uI[0]; PI2[m+2] =   uI[m];
    PI3[0] =   pI[0]; PI3[m+2] =   pI[m];
  }
  else
  {
    W1[0]   = rho[m-3];W1[1]   = rho[m-2];W1[2]   = rho[m-1];
    W2[0]   = mom[m-3];W2[1]   = mom[m-2];W2[2]   = mom[m-1];
    W3[0]   = ene[m-3];W3[1]   = ene[m-2];W3[2]   = ene[m-1];
    W1[m+3] = rho[0];  W1[m+4] = rho[1];  W1[m+5] = rho[2];
    W2[m+3] = mom[0];  W2[m+4] = mom[1];  W2[m+5] = mom[2];
    W3[m+3] = ene[0];  W3[m+4] = ene[1];  W3[m+5] = ene[2];
    PI1[0] = rhoI[m-1]; PI1[m+2] = rhoI[1];
    PI2[0] =   uI[m-1]; PI2[m+2] =   uI[1];
    PI3[0] =   pI[m-1]; PI3[m+2] =   pI[1];
  }

  for(j = 0; j < m+6; ++j)
  {
    u[j] = W2[j]/W1[j];
    p[j] = (W3[j] - 0.5*W2[j]*u[j])*(gamma-1.0);
    H[j] = 0.5*u[j]*u[j] + gamma*p[j]/W1[j]/(gamma-1.0);
  }

  for(j = 3; j < m+4; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 6; ++k)
    {
      Q1[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-3+k]-0.5*W1[j-3+k]*u_star)*u_star - W3[j-3+k]) / c_star/c_star + W1[j-3+k];
    }

    if(WENOD)
    {
      local_WENO_5_inter_d(h, Q1, DQQ[0]);
      local_WENO_5_inter_d(h, Q2, DQQ[1]);
      local_WENO_5_inter_d(h, Q3, DQQ[2]);
    }
    else
    {
      local_WENO_5_inter(h, Q1);
      local_WENO_5_inter(h, Q2);
      local_WENO_5_inter(h, Q3);
      DQQ[0][0] = (15*(Q1[3] - Q1[2]) + Q1[1] - Q1[4]) / (12*h);
      DQQ[1][0] = (15*(Q2[3] - Q2[2]) + Q2[1] - Q2[4]) / (12*h);
      DQQ[2][0] = (15*(Q3[3] - Q3[2]) + Q3[1] - Q3[4]) / (12*h);
      DQQ[0][1] = DQQ[0][0];
      DQQ[1][1] = DQQ[1][0];
      DQQ[2][1] = DQQ[2][0];
    }

  //=====Recomposition========
    rho_L[j-3] = Q1[6] + Q2[6] + Q3[6];
    rho_R[j-3] = Q1[7] + Q2[7] + Q3[7];
    D_rho_L[j-3] = DQQ[0][0] + DQQ[1][0] + DQQ[2][0];
    D_rho_R[j-3] = DQQ[0][1] + DQQ[1][1] + DQQ[2][1];

    u_L[j-3] = u_star*rho_L[j-3] + c_star*(Q3[6] - Q1[6]);
    u_R[j-3] = u_star*rho_R[j-3] + c_star*(Q3[7] - Q1[7]);
    D_u_L[j-3] = u_star*D_rho_L[j-3] + c_star*(DQQ[2][0] - DQQ[0][0]);
    D_u_R[j-3] = u_star*D_rho_R[j-3] + c_star*(DQQ[2][1] - DQQ[0][1]);

    p_L[j-3] = H_star*(Q1[6]+Q3[6]) + u_star*c_star*(Q3[6]-Q1[6]) + 0.5*u_star*u_star*Q2[6];
    p_R[j-3] = H_star*(Q1[7]+Q3[7]) + u_star*c_star*(Q3[7]-Q1[7]) + 0.5*u_star*u_star*Q2[7];
    D_p_L[j-3] = H_star*(DQQ[0][0]+DQQ[2][0]) + u_star*c_star*(DQQ[2][0]-DQQ[0][0]) + 0.5*u_star*u_star*DQQ[1][0];
    D_p_R[j-3] = H_star*(DQQ[0][1]+DQQ[2][1]) + u_star*c_star*(DQQ[2][1]-DQQ[0][1]) + 0.5*u_star*u_star*DQQ[1][1];


    u_L[j-3] = u_L[j-3] / rho_L[j-3];
    u_R[j-3] = u_R[j-3] / rho_R[j-3];
    p_L[j-3] = (p_L[j-3] - 0.5*u_L[j-3]*u_L[j-3]*rho_L[j-3])*(gamma-1.0);
    p_R[j-3] = (p_R[j-3] - 0.5*u_R[j-3]*u_R[j-3]*rho_R[j-3])*(gamma-1.0);

    D_u_L[j-3] = (D_u_L[j-3] - u_L[j-3]*D_rho_L[j-3]) / rho_L[j-3];
    D_u_R[j-3] = (D_u_R[j-3] - u_R[j-3]*D_rho_R[j-3]) / rho_R[j-3];
    D_p_L[j-3] = (D_p_L[j-3] - 0.5*D_rho_L[j-3]*u_L[j-3]*u_L[j-3] - rho_L[j-3]*u_L[j-3]*D_u_L[j-3])*(gamma-1.0);
    D_p_R[j-3] = (D_p_R[j-3] - 0.5*D_rho_R[j-3]*u_R[j-3]*u_R[j-3] - rho_R[j-3]*u_R[j-3]*D_u_R[j-3])*(gamma-1.0);


    flag = 0;

    deltaP = threshold*fabs(rho_L[j-3] - rho_R[j-3]);
    DP = 1e-5+fabs(10*(W1[j]-W1[j-1]) + 5*(W1[j-2]-W1[j+1]) + (W1[j+2]-W1[j-3]));
    if(deltaP > DP) flag = 4*limiter;

    for(k = 0; k < 6; ++k)
      Q1[k] = Entp[j-3+k];
    local_WENO_5_inter(h, Q3);
    deltaP = threshold*fabs(Q1[6] - Q1[7]);
    DP = 1e-5+fabs(10*(Q1[3]-Q1[2]) + 5*(Q1[1]-Q1[4]) + (Q1[5]-Q1[0]));
    if(deltaP > DP) flag = 5*limiter;

    if(j-3)
      trouble[j-4] += flag;
    if(m+3-j)
      trouble[j-3] += flag;
  }


  for(j = 0; j < m; ++j)
  {
    if(trouble[j])
    {
      Qn1 = W1[j+2];
      Q0  = W1[j+3];
      Qp1 = W1[j+4];

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

      Qn1 = u[j+2];
      Q0  = u[j+3];
      Qp1 = u[j+4];
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

      Qn1 = p[j+2];
      Q0  = p[j+3];
      Qp1 = p[j+4];
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







void WENO_50
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  double const threshold = running_info[6];

  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double Q1[8], Q2[8], Q3[8];
  double SL, SR, S0, Qn2, Qn1, Q0, Qp1, Qp2, DQQ[3][2];
  double H[m+6], u[m+6], p[m+6];
  double deltaP, DP, flag, trouble[m];

  double H_star, u_star, c_star;
  for(j = 0; j < m; ++j)
  {
    W1[j+3] = rho[j];
    W2[j+3] = mom[j];
    W3[j+3] = ene[j];
  }

  if(bod < 0)
  {
    W1[0]   = rho[2];  W1[1]   = rho[1];  W1[2]   = rho[0];
    W2[0]   =-mom[2];  W2[1]   =-mom[1];  W2[2]   =-mom[0];
    W3[0]   = ene[2];  W3[1]   = ene[1];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-2];W1[m+5] = rho[m-3];
    W2[m+3] =-mom[m-1];W2[m+4] =-mom[m-2];W2[m+5] =-mom[m-3];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-2];W3[m+5] = ene[m-3];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];  W1[2]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];  W2[2]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-1];W1[m+5] = rho[m-1];
    W2[m+3] = mom[m-1];W2[m+4] = mom[m-1];W2[m+5] = mom[m-1];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-1];W3[m+5] = ene[m-1];
  }
  else
  {
    W1[0]   = rho[m-3];W1[1]   = rho[m-2];W1[2]   = rho[m-1];
    W2[0]   = mom[m-3];W2[1]   = mom[m-2];W2[2]   = mom[m-1];
    W3[0]   = ene[m-3];W3[1]   = ene[m-2];W3[2]   = ene[m-1];
    W1[m+3] = rho[0];  W1[m+4] = rho[1];  W1[m+5] = rho[2];
    W2[m+3] = mom[0];  W2[m+4] = mom[1];  W2[m+5] = mom[2];
    W3[m+3] = ene[0];  W3[m+4] = ene[1];  W3[m+5] = ene[2];
  }

  for(j = 0; j < m+6; ++j)
  {
    u[j] = W2[j]/W1[j];
    p[j] = (W3[j] - 0.5*W2[j]*u[j])*(gamma-1.0);
    H[j] = 0.5*u[j]*u[j] + gamma*p[j]/W1[j]/(gamma-1.0);
  }

  for(j = 3; j < m+4; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 6; ++k)
    {
      Q1[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-3+k]-0.5*W1[j-3+k]*u_star)*u_star - W3[j-3+k]) / c_star/c_star + W1[j-3+k];
    }

    if(WENOD)
    {
      local_WENO_5_inter_d(h, Q1, DQQ[0]);
      local_WENO_5_inter_d(h, Q2, DQQ[1]);
      local_WENO_5_inter_d(h, Q3, DQQ[2]);
    }
    else
    {
      local_WENO_5_inter(h, Q1);
      local_WENO_5_inter(h, Q2);
      local_WENO_5_inter(h, Q3);
      DQQ[0][0] = (15*(Q1[3] - Q1[2]) + Q1[1] - Q1[4]) / (12*h);
      DQQ[1][0] = (15*(Q2[3] - Q2[2]) + Q2[1] - Q2[4]) / (12*h);
      DQQ[2][0] = (15*(Q3[3] - Q3[2]) + Q3[1] - Q3[4]) / (12*h);
      DQQ[0][1] = DQQ[0][0];
      DQQ[1][1] = DQQ[1][0];
      DQQ[2][1] = DQQ[2][0];
    }


  //=====Recomposition========
    rho_L[j-3] = Q1[6] + Q2[6] + Q3[6];
    rho_R[j-3] = Q1[7] + Q2[7] + Q3[7];
    D_rho_L[j-3] = DQQ[0][0] + DQQ[1][0] + DQQ[2][0];
    D_rho_R[j-3] = DQQ[0][1] + DQQ[1][1] + DQQ[2][1];

    u_L[j-3] = u_star*rho_L[j-3] + c_star*(Q3[6] - Q1[6]);
    u_R[j-3] = u_star*rho_R[j-3] + c_star*(Q3[7] - Q1[7]);
    D_u_L[j-3] = u_star*D_rho_L[j-3] + c_star*(DQQ[2][0] - DQQ[0][0]);
    D_u_R[j-3] = u_star*D_rho_R[j-3] + c_star*(DQQ[2][1] - DQQ[0][1]);

    p_L[j-3] = H_star*(Q1[6]+Q3[6]) + u_star*c_star*(Q3[6]-Q1[6]) + 0.5*u_star*u_star*Q2[6];
    p_R[j-3] = H_star*(Q1[7]+Q3[7]) + u_star*c_star*(Q3[7]-Q1[7]) + 0.5*u_star*u_star*Q2[7];
    D_p_L[j-3] = H_star*(DQQ[0][0]+DQQ[2][0]) + u_star*c_star*(DQQ[2][0]-DQQ[0][0]) + 0.5*u_star*u_star*DQQ[1][0];
    D_p_R[j-3] = H_star*(DQQ[0][1]+DQQ[2][1]) + u_star*c_star*(DQQ[2][1]-DQQ[0][1]) + 0.5*u_star*u_star*DQQ[1][1];


    u_L[j-3] = u_L[j-3] / rho_L[j-3];
    u_R[j-3] = u_R[j-3] / rho_R[j-3];
    p_L[j-3] = (p_L[j-3] - 0.5*u_L[j-3]*u_L[j-3]*rho_L[j-3])*(gamma-1.0);
    p_R[j-3] = (p_R[j-3] - 0.5*u_R[j-3]*u_R[j-3]*rho_R[j-3])*(gamma-1.0);

    D_u_L[j-3] = (D_u_L[j-3] - u_L[j-3]*D_rho_L[j-3]) / rho_L[j-3];
    D_u_R[j-3] = (D_u_R[j-3] - u_R[j-3]*D_rho_R[j-3]) / rho_R[j-3];
    D_p_L[j-3] = (D_p_L[j-3] - 0.5*D_rho_L[j-3]*u_L[j-3]*u_L[j-3] - rho_L[j-3]*u_L[j-3]*D_u_L[j-3])*(gamma-1.0);
    D_p_R[j-3] = (D_p_R[j-3] - 0.5*D_rho_R[j-3]*u_R[j-3]*u_R[j-3] - rho_R[j-3]*u_R[j-3]*D_u_R[j-3])*(gamma-1.0);



    flag = 0;

    deltaP = threshold*fabs(rho_L[j-3] - rho_R[j-3]);
    DP = 1e-5+fabs(10*(W1[j]-W1[j-1]) + 5*(W1[j-2]-W1[j+1]) + (W1[j+2]-W1[j-3]));
    if(deltaP > DP) flag = 4*limiter;

    if(j-3)
      trouble[j-4] += flag;
    if(m+3-j)
      trouble[j-3] += flag;
  }


  for(j = 0; j < m; ++j)
  {
    if(trouble[j])
    {
      Qn1 = W1[j+2];
      Q0  = W1[j+3];
      Qp1 = W1[j+4];

      SL = Q0 - Qn1; SR = Qp1 - Q0;
      if(SL*SR < 0.0)
      {
	D_rho_R[j] = 0.0;
	rho_R[j] = Q0;
	rho_L[j+1] = Q0;
	D_rho_L[j+1] = D_rho_R[j];
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	rho_R[j] = Q0 - 0.5*S0;
	D_rho_R[j] = S0 / h;
	rho_L[j+1] = Q0 + 0.5*S0;
	D_rho_L[j+1] = D_rho_R[j];
      }

      Qn1 = u[j+2];
      Q0  = u[j+3];
      Qp1 = u[j+4];
      SL = Q0 - Qn1; SR = Qp1 - Q0;
      if(SL*SR < 0.0)
      {
	D_u_R[j] = 0.0;
	u_R[j] = Q0;
	u_L[j+1] = Q0;
	D_u_L[j+1] = D_u_R[j];
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	u_R[j] = Q0 - 0.5*S0;
	D_u_R[j] = S0 / h;
	u_L[j+1] = Q0 + 0.5*S0;
	D_u_L[j+1] = D_u_R[j];
      }

      Qn1 = p[j+2];
      Q0  = p[j+3];
      Qp1 = p[j+4];
      SL = Q0 - Qn1; SR = Qp1 - Q0;
      if(SL*SR < 0.0)
      {
	D_p_R[j] = 0.0;
	p_R[j] = Q0;
	p_L[j+1] = Q0;
	D_p_L[j+1] = D_p_R[j];
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	p_R[j] = Q0 - 0.5*S0;
	D_p_R[j] = S0 / h;
	p_L[j+1] = Q0 + 0.5*S0;
	D_p_L[j+1] = D_p_R[j];
      }
    }
  }
}




void WENO_5
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  double const threshold = running_info[6];

  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double Q1[8], Q2[8], Q3[8];
  double SL, SR, S0, Qn2, Qn1, Qp1, Qp2, DQQ[3][2];
  double H[m+6], u, p;

  double H_star, u_star, c_star;
  for(j = 0; j < m; ++j)
  {
    W1[j+3] = rho[j];
    W2[j+3] = mom[j];
    W3[j+3] = ene[j];
  }

  if(bod < 0)
  {
    W1[0]   = rho[2];  W1[1]   = rho[1];  W1[2]   = rho[0];
    W2[0]   =-mom[2];  W2[1]   =-mom[1];  W2[2]   =-mom[0];
    W3[0]   = ene[2];  W3[1]   = ene[1];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-2];W1[m+5] = rho[m-3];
    W2[m+3] =-mom[m-1];W2[m+4] =-mom[m-2];W2[m+5] =-mom[m-3];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-2];W3[m+5] = ene[m-3];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];  W1[2]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];  W2[2]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-1];W1[m+5] = rho[m-1];
    W2[m+3] = mom[m-1];W2[m+4] = mom[m-1];W2[m+5] = mom[m-1];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-1];W3[m+5] = ene[m-1];
  }
  else
  {
    W1[0]   = rho[m-3];W1[1]   = rho[m-2];W1[2]   = rho[m-1];
    W2[0]   = mom[m-3];W2[1]   = mom[m-2];W2[2]   = mom[m-1];
    W3[0]   = ene[m-3];W3[1]   = ene[m-2];W3[2]   = ene[m-1];
    W1[m+3] = rho[0];  W1[m+4] = rho[1];  W1[m+5] = rho[2];
    W2[m+3] = mom[0];  W2[m+4] = mom[1];  W2[m+5] = mom[2];
    W3[m+3] = ene[0];  W3[m+4] = ene[1];  W3[m+5] = ene[2];
  }
  /*
  W1[m+3] = 0.2*(cos(25.0+5.0*h) - cos(25.0));
  W1[m+3] = 1.0 + 0.2*W1[m+3];
  W1[m+4] = 0.2*(cos(25.0+10.0*h) - cos(25.0+5.0*h));
  W1[m+4] = 1.0 + 0.2*W1[m+4];
  W1[m+5] = 0.2*(cos(25.0+15.0*h) - cos(25.0+10.0*h));
  W1[m+5] = 1.0 + 0.2*W1[m+5];
  */

  for(j = 0; j < m+6; ++j)
  {
    u = W2[j]/W1[j];
    p = (W3[j] - 0.5*W2[j]*u)*(gamma-1.0);
    H[j] = 0.5*u*u + gamma*p/W1[j]/(gamma-1.0);
  }

  for(j = 3; j < m+4; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 6; ++k)
    {
      Q1[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-3+k]-0.5*W1[j-3+k]*u_star)*u_star - W3[j-3+k]) / c_star/c_star + W1[j-3+k];
    }

    if(WENOD)
    {
      local_WENO_5_inter_d(h, Q1, DQQ[0]);
      local_WENO_5_inter_d(h, Q2, DQQ[1]);
      local_WENO_5_inter_d(h, Q3, DQQ[2]);
    }
    else
    {
      local_WENO_5_inter(h, Q1);
      local_WENO_5_inter(h, Q2);
      local_WENO_5_inter(h, Q3);
      DQQ[0][0] = (15*(Q1[3] - Q1[2]) + Q1[1] - Q1[4]) / (12*h);
      DQQ[1][0] = (15*(Q2[3] - Q2[2]) + Q2[1] - Q2[4]) / (12*h);
      DQQ[2][0] = (15*(Q3[3] - Q3[2]) + Q3[1] - Q3[4]) / (12*h);
      DQQ[0][1] = DQQ[0][0];
      DQQ[1][1] = DQQ[1][0];
      DQQ[2][1] = DQQ[2][0];
    }


  //=====Recomposition========
    rho_L[j-3] = Q1[6] + Q2[6] + Q3[6];
    rho_R[j-3] = Q1[7] + Q2[7] + Q3[7];
    D_rho_L[j-3] = DQQ[0][0] + DQQ[1][0] + DQQ[2][0];
    D_rho_R[j-3] = DQQ[0][1] + DQQ[1][1] + DQQ[2][1];

    u_L[j-3] = u_star*rho_L[j-3] + c_star*(Q3[6] - Q1[6]);
    u_R[j-3] = u_star*rho_R[j-3] + c_star*(Q3[7] - Q1[7]);
    D_u_L[j-3] = u_star*D_rho_L[j-3] + c_star*(DQQ[2][0] - DQQ[0][0]);
    D_u_R[j-3] = u_star*D_rho_R[j-3] + c_star*(DQQ[2][1] - DQQ[0][1]);

    p_L[j-3] = H_star*(Q1[6]+Q3[6]) + u_star*c_star*(Q3[6]-Q1[6]) + 0.5*u_star*u_star*Q2[6];
    p_R[j-3] = H_star*(Q1[7]+Q3[7]) + u_star*c_star*(Q3[7]-Q1[7]) + 0.5*u_star*u_star*Q2[7];
    D_p_L[j-3] = H_star*(DQQ[0][0]+DQQ[2][0]) + u_star*c_star*(DQQ[2][0]-DQQ[0][0]) + 0.5*u_star*u_star*DQQ[1][0];
    D_p_R[j-3] = H_star*(DQQ[0][1]+DQQ[2][1]) + u_star*c_star*(DQQ[2][1]-DQQ[0][1]) + 0.5*u_star*u_star*DQQ[1][1];


    u_L[j-3] = u_L[j-3] / rho_L[j-3];
    u_R[j-3] = u_R[j-3] / rho_R[j-3];
    p_L[j-3] = (p_L[j-3] - 0.5*u_L[j-3]*u_L[j-3]*rho_L[j-3])*(gamma-1.0);
    p_R[j-3] = (p_R[j-3] - 0.5*u_R[j-3]*u_R[j-3]*rho_R[j-3])*(gamma-1.0);

    D_u_L[j-3] = (D_u_L[j-3] - u_L[j-3]*D_rho_L[j-3]) / rho_L[j-3];
    D_u_R[j-3] = (D_u_R[j-3] - u_R[j-3]*D_rho_R[j-3]) / rho_R[j-3];
    D_p_L[j-3] = (D_p_L[j-3] - 0.5*D_rho_L[j-3]*u_L[j-3]*u_L[j-3] - rho_L[j-3]*u_L[j-3]*D_u_L[j-3])*(gamma-1.0);
    D_p_R[j-3] = (D_p_R[j-3] - 0.5*D_rho_R[j-3]*u_R[j-3]*u_R[j-3] - rho_R[j-3]*u_R[j-3]*D_u_R[j-3])*(gamma-1.0);
  }
}








void WENO_5_noD
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[])
{
  int j, k;
  int const bod = (int)running_info[3];
  int const WENOD = (int)running_info[4];
  int const limiter = (int)running_info[5];
  double const threshold = running_info[6];


  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double PI1[m+3], PI2[m+3], PI3[m+3];
  double Q1[8], Q2[8], Q3[8];
  double u[m+6], p[m+6], H[m+6], Entp[m+6];
  double EntpL, EntpR;
  double Qn1, Q0, Qp1, SL, SR, S0;
  double deltaP, DP;
  int flag, trouble[m];

  double H_star, u_star, c_star;
  for(j = 0; j < m; ++j)
  {
    W1[j+3] = rho[j];
    W2[j+3] = mom[j];
    W3[j+3] = ene[j];
    trouble[j] = 0;
  }

  if(bod < 0)
  {
    W1[0]   = rho[2];  W1[1]   = rho[1];  W1[2]   = rho[0];
    W2[0]   =-mom[2];  W2[1]   =-mom[1];  W2[2]   =-mom[0];
    W3[0]   = ene[2];  W3[1]   = ene[1];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-2];W1[m+5] = rho[m-3];
    W2[m+3] =-mom[m-1];W2[m+4] =-mom[m-2];W2[m+5] =-mom[m-3];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-2];W3[m+5] = ene[m-3];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];  W1[1]   = rho[0];  W1[2]   = rho[0];
    W2[0]   = mom[0];  W2[1]   = mom[0];  W2[2]   = mom[0];
    W3[0]   = ene[0];  W3[1]   = ene[0];  W3[2]   = ene[0];
    W1[m+3] = rho[m-1];W1[m+4] = rho[m-1];W1[m+5] = rho[m-1];
    W2[m+3] = mom[m-1];W2[m+4] = mom[m-1];W2[m+5] = mom[m-1];
    W3[m+3] = ene[m-1];W3[m+4] = ene[m-1];W3[m+5] = ene[m-1];
  }
  else
  {
    W1[0]   = rho[m-3];W1[1]   = rho[m-2];W1[2]   = rho[m-1];
    W2[0]   = mom[m-3];W2[1]   = mom[m-2];W2[2]   = mom[m-1];
    W3[0]   = ene[m-3];W3[1]   = ene[m-2];W3[2]   = ene[m-1];
    W1[m+3] = rho[0];  W1[m+4] = rho[1];  W1[m+5] = rho[2];
    W2[m+3] = mom[0];  W2[m+4] = mom[1];  W2[m+5] = mom[2];
    W3[m+3] = ene[0];  W3[m+4] = ene[1];  W3[m+5] = ene[2];
  }

  for(j = 0; j < m+6; ++j)
  {
    u[j] = W2[j]/W1[j];
    p[j] = (W3[j] - 0.5*W2[j]*u[j])*(gamma-1.0);
    H[j] = 0.5*u[j]*u[j] + gamma*p[j]/W1[j]/(gamma-1.0);
  }

  for(j = 3; j < m+4; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_star = sqrt((gamma-1.0)*(H_star - 0.5*u_star*u_star));

  //=========Charactoristic Decomposition=========
    for(k = 0; k < 6; ++k)
    {
      Q1[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = W1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= W2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= W3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((W2[j-3+k]-0.5*W1[j-3+k]*u_star)*u_star - W3[j-3+k]) / c_star/c_star + W1[j-3+k];
    }

    local_WENO_5_inter(h, Q1);
    local_WENO_5_inter(h, Q2);
    local_WENO_5_inter(h, Q3);

  //=====Recomposition========
    rho_L[j-3] = Q1[6] + Q2[6] + Q3[6];
    rho_R[j-3] = Q1[7] + Q2[7] + Q3[7];

    u_L[j-3] = u_star*rho_L[j-3] + c_star*(Q3[6] - Q1[6]);
    u_R[j-3] = u_star*rho_R[j-3] + c_star*(Q3[7] - Q1[7]);

    p_L[j-3] = H_star*(Q1[6]+Q3[6]) + u_star*c_star*(Q3[6]-Q1[6]) + 0.5*u_star*u_star*Q2[6];
    p_R[j-3] = H_star*(Q1[7]+Q3[7]) + u_star*c_star*(Q3[7]-Q1[7]) + 0.5*u_star*u_star*Q2[7];


    u_L[j-3] = u_L[j-3] / rho_L[j-3];
    u_R[j-3] = u_R[j-3] / rho_R[j-3];
    p_L[j-3] = (p_L[j-3] - 0.5*u_L[j-3]*u_L[j-3]*rho_L[j-3])*(gamma-1.0);
    p_R[j-3] = (p_R[j-3] - 0.5*u_R[j-3]*u_R[j-3]*rho_R[j-3])*(gamma-1.0);


    flag = 0;

    deltaP = threshold*fabs(rho_L[j-3] - rho_R[j-3]);
    DP = 1e-5+fabs(10*(W1[j]-W1[j-1]) + 5*(W1[j-2]-W1[j+1]) + (W1[j+2]-W1[j-3]));
    if(deltaP > DP) flag = 4*limiter;

    for(k = 0; k < 6; ++k)
      Q1[k] = Entp[j-3+k];
    local_WENO_5_inter(h, Q3);
    deltaP = threshold*fabs(Q1[6] - Q1[7]);
    DP = 1e-5+fabs(10*(Q1[3]-Q1[2]) + 5*(Q1[1]-Q1[4]) + (Q1[5]-Q1[0]));
    if(deltaP > DP) flag = 5*limiter;

    if(j-3)
      trouble[j-3] += flag;
    if(m+3-j)
      trouble[j-2] += flag;
  }


  for(j = 0; j < m; ++j)
  {
    if(trouble[j])
    {
      Qn1 = W1[j+2];
      Q0  = W1[j+3];
      Qp1 = W1[j+4];

      SL = Q0 - Qn1; SR = Qp1 - Q0;;
      if(SL*SR < 0.0)
      {
	rho_R[j] = Q0;
	rho_L[j+1] = Q0;
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	rho_R[j] = Q0 - 0.5*S0;
	rho_L[j+1] = Q0 + 0.5*S0;
      }

      Qn1 = u[j+2];
      Q0  = u[j+3];
      Qp1 = u[j+4];
      SL = Q0 - Qn1; SR = Qp1 - Q0;;
      if(SL*SR < 0.0)
      {
	u_R[j] = Q0;
	u_L[j+1] = Q0;
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	u_R[j] = Q0 - 0.5*S0;
	u_L[j+1] = Q0 + 0.5*S0;
      }

      Qn1 = p[j+2];
      Q0  = p[j+3];
      Qp1 = p[j+4];
      SL = Q0 - Qn1; SR = Qp1 - Q0;
      if(SL*SR < 0.0)
      {
	p_R[j] = Q0;
	p_L[j+1] = Q0;
      }
      else
      {
	S0 = alp2*((fabs(SL) < fabs(SR))? SL : SR);
	p_R[j] = Q0 - 0.5*S0;
	p_L[j+1] = Q0 + 0.5*S0;
      }
    }
  }
}
