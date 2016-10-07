#include <math.h>

#include "reconstruction.h"

#include "de_re_composition.h"


void WENO_30
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];
  int j, k;
  int l = 4, half_l = 2;

  double W1[m+l], W2[m+l], W3[m+l]; //W1=rho, W2=mom, W3=ene
  double Q1[l+2], Q2[l+2], Q3[l+2];
  double SL, SR, S0, Qn2, Qn1, Q0, Qp1, Qp2, DQQ[3][2];
  double QL[3], QR[3], PL[3], PR[3];
  double H[m+l], u[m+l], p[m+l];
  double deltaP, DP, flag, trouble[m];

  double H_star, u_star, c_square, c_star, gamma1 = gamma-1.0;
  for(j = 0; j < m; ++j)
  {
    W1[j+half_l] = rho[j];
    W2[j+half_l] = mom[j];
    W3[j+half_l] = ene[j];
  }

  if(bod < 0)
  {
    W1[0]   = rho[1];        W1[1]   = rho[0];
    W2[0]   =-mom[1];        W2[1]   =-mom[0];
    W3[0]   = ene[1];        W3[1]   = ene[0];
    W1[m+half_l] = rho[m-1]; W1[m+half_l+1] = rho[m-2];
    W2[m+half_l] =-mom[m-1]; W2[m+half_l+1] =-mom[m-2];
    W3[m+half_l] = ene[m-1]; W3[m+half_l+1] = ene[m-2];
  }
  else if(bod > 0)
  {
    W1[0]   = rho[0];        W1[1]   = rho[0];
    W2[0]   = mom[0];        W2[1]   = mom[0];
    W3[0]   = ene[0];        W3[1]   = ene[0];
    W1[m+half_l] = rho[m-1]; W1[m+half_l+1] = rho[m-1];
    W2[m+half_l] = mom[m-1]; W2[m+half_l+1] = mom[m-1];
    W3[m+half_l] = ene[m-1]; W3[m+half_l+1] = ene[m-1];
  }
  else
  {
    W1[0]   = rho[m-2];     W1[1]   = rho[m-1];
    W2[0]   = mom[m-2];     W2[1]   = mom[m-1];
    W3[0]   = ene[m-2];     W3[1]   = ene[m-1];
    W1[m+half_l] = rho[0];  W1[m+half_l+1] = rho[1];
    W2[m+half_l] = mom[0];  W2[m+half_l+1] = mom[1];
    W3[m+half_l] = ene[0];  W3[m+half_l+1] = ene[1];
  }

  for(j = 0; j < m+l; ++j)
  {
    u[j] = W2[j]/W1[j];
    p[j] = (W3[j] - 0.5*W2[j]*u[j])*(gamma-1.0);
    H[j] = 0.5*u[j]*u[j] + gamma*p[j]/W1[j]/(gamma-1.0);
  }

  for(j = half_l; j < m+l; ++j)
  {
    u_star = (W2[j]/sqrt(W1[j]) + W2[j-1]/sqrt(W1[j-1]))/(sqrt(W1[j]) + sqrt(W1[j-1]));
    H_star = (sqrt(W1[j])* H[j] + sqrt(W1[j-1])* H[j-1])/(sqrt(W1[j]) + sqrt(W1[j-1]));
    c_square = gamma1*(H_star - 0.5*u_star*u_star);
    c_star = sqrt(c_square);

  //=========Charactoristic Decomposition=========
    decomposition(H_star, u_star, c_square, c_star, gamma1, m, j, l, half_l, decomp, W1, W2, W3, Q1, Q2, Q3);

    local_WENO_3_inter_Z(h, Q1);
    local_WENO_3_inter_Z(h, Q2);
    local_WENO_3_inter_Z(h, Q3);
    DQQ[0][0] = (Q1[2] - Q1[1]) / h;
    DQQ[1][0] = (Q2[2] - Q2[1]) / h;
    DQQ[2][0] = (Q3[2] - Q3[1]) / h;
    DQQ[0][1] = DQQ[0][0];
    DQQ[1][1] = DQQ[1][0];
    DQQ[2][1] = DQQ[2][0];


  //=====Recomposition========
    QL[0] = Q1[l]; QR[0] = Q1[l+1];
    QL[1] = Q2[l]; QR[1] = Q2[l+1];
    QL[2] = Q3[l]; QR[2] = Q3[l+1];
    recomposition(H_star, u_star, c_star, decomp, QL, QR, PL, PR);
    rho_L[j-half_l] = PL[0]; rho_R[j-half_l] = PR[0];
    u_L[j-half_l] = PL[1];   u_R[j-half_l] = PR[1];
    p_L[j-half_l] = PL[2];   p_R[j-half_l] = PR[2];
    QL[0] = DQQ[0][0]; QR[0] = DQQ[0][1];
    QL[1] = DQQ[1][0]; QR[1] = DQQ[1][1];
    QL[2] = DQQ[2][0]; QR[2] = DQQ[2][1];
    recomposition(H_star, u_star, c_star, decomp, QL, QR, PL, PR);
    D_rho_L[j-half_l] = PL[0]; D_rho_R[j-half_l] = PR[0];
    D_u_L[j-half_l] = PL[1];   D_u_R[j-half_l] = PR[1];
    D_p_L[j-half_l] = PL[2];   D_p_R[j-half_l] = PR[2];

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
