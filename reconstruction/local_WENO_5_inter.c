#include <math.h>


void local_WENO_5_inter_d(double h, double Q[8], double DQQ[2])
{
  double c[3][3], d[3][3], e[3][3], f[3][3], eps = 1e-6;
  c[0][0] = 1.0/3.0;
  c[0][1] = 5.0/6.0;
  c[0][2] = -1.0/6.0;
  c[1][0] = c[0][2];
  c[1][1] = c[0][1];
  c[1][2] = c[0][0];
  c[2][0] = 1.0/3.0;
  c[2][1] = -7.0/6.0;
  c[2][2] = 11.0/6.0;
  d[1][0] = c[0][0];
  d[1][1] = c[0][1];
  d[1][2] = c[0][2];
  d[2][0] = c[1][0];
  d[2][1] = c[1][1];
  d[2][2] = c[1][2];
  d[0][0] = 11.0/6.0;
  d[0][1] = -7.0/6.0;
  d[0][2] = 1.0/3.0;
  e[0][0] = -1.0;
  e[0][1] = 1.0;
  e[0][2] = 0.0;
  e[1][0] = -e[0][2];
  e[1][1] = -e[0][1];
  e[1][2] = -e[0][0];
  e[2][0] = 1.0;
  e[2][1] = -3.0;
  e[2][2] = 2.0;
  f[1][0] = e[0][0];
  f[1][1] = e[0][1];
  f[1][2] = e[0][2];
  f[2][0] = e[1][0];
  f[2][1] = e[1][1];
  f[2][2] = e[1][2];
  f[0][0] = -e[2][2];
  f[0][1] = -e[2][1];
  f[0][2] = -e[2][0];

  double Dneg[3] = {0.3, 0.6, 0.1};//dR
  double Dpos[3] = {0.1, 0.6, 0.3};//dL
  double dDneg[3] = {0.5, 5.0/12.0, 1.0/12.0};//dDR
  double dDpos[3] = {1.0/12.0, 5.0/12.0, 0.5};//dDL

  double alp[3], omg[3], beta[3], sum;
  double Qneg[3], Qpos[3], dQneg[3], dQpos[3];


  //=====Reconstruction for The Negative Limit======
    Qneg[2] = c[2][0]*Q[0] + c[2][1]*Q[1] + c[2][2]*Q[2];
    Qneg[1] = c[1][0]*Q[1] + c[1][1]*Q[2] + c[1][2]*Q[3];
    Qneg[0] = c[0][0]*Q[2] + c[0][1]*Q[3] + c[0][2]*Q[4];

    dQneg[2] = (e[2][0]*Q[0] + e[2][1]*Q[1] + e[2][2]*Q[2])/h;
    dQneg[1] = (               e[1][1]*Q[2] + e[1][2]*Q[3])/h;
    dQneg[0] = dQneg[1];
    //dQneg[0] = (e[0][0]*Q1[2] + e[0][1]*Q1[3]                )/h;

    beta[0] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[1] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    beta[2] = (Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
    alp[0] = dDneg[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = dDneg[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = dDneg[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    DQQ[0] = dQneg[0]*omg[0] + dQneg[1]*omg[1] + dQneg[2]*omg[2];

    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[2] - 4.0*Q[3]     + Q[4])*(3.0*Q[2] - 4.0*Q[3]     + Q[4]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[1]                - Q[3])*(    Q[1]                - Q[3]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[0] - 4.0*Q[1] + 3.0*Q[2])*(    Q[0] - 4.0*Q[1] + 3.0*Q[2]);
    alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dneg[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[6] = Qneg[0]*omg[0] + Qneg[1]*omg[1] + Qneg[2]*omg[2];


  //=====Reconstruction for The Positive Limit======
    Qpos[2] = d[2][0]*Q[1] + d[2][1]*Q[2] + d[2][2]*Q[3];
    Qpos[1] = d[1][0]*Q[2] + d[1][1]*Q[3] + d[1][2]*Q[4];
    Qpos[0] = d[0][0]*Q[3] + d[0][1]*Q[4] + d[0][2]*Q[5];

    dQpos[0] = (f[0][0]*Q[3] + f[0][1]*Q[4] + f[0][2]*Q[5])/h;
    dQpos[1] = (f[1][0]*Q[2] + f[1][1]*Q[3]               )/h;
    dQpos[2] = dQpos[1];
    //dQpos[0] = (f[0][0]*Q1[3] + f[0][1]*Q1[4] + f[0][2]*Q1[5])/h;

    beta[0] = (Q[3] - 2.0*Q[4] + Q[5])*(Q[3] - 2.0*Q[4] + Q[5]);
    beta[1] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[2] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    alp[0] = dDpos[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = dDpos[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = dDpos[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    DQQ[1] = dQpos[0]*omg[0] + dQpos[1]*omg[1] + dQpos[2]*omg[2];

    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[3] - 4.0*Q[4]     + Q[5])*(3.0*Q[3] - 4.0*Q[4]     + Q[5]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[2]                - Q[4])*(    Q[2]                - Q[4]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[1] - 4.0*Q[2] + 3.0*Q[3])*(    Q[1] - 4.0*Q[2] + 3.0*Q[3]);
    alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dpos[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[7] = Qpos[0]*omg[0] + Qpos[1]*omg[1] + Qpos[2]*omg[2];
}




void local_WENO_5_interleft(double h, double Q[8])
{
  double c[3][3], eps = 1e-6;
  c[0][0] = 1.0/3.0;
  c[0][1] = 5.0/6.0;
  c[0][2] = -1.0/6.0;
  c[1][0] = c[0][2];
  c[1][1] = c[0][1];
  c[1][2] = c[0][0];
  c[2][0] = 1.0/3.0;
  c[2][1] = -7.0/6.0;
  c[2][2] = 11.0/6.0;

  double Dneg[3] = {0.3, 0.6, 0.1};//dR

  double alp[3], omg[3], beta[3], sum;
  double Qneg[3], dQneg[3];


  //=====Reconstruction for The Negative Limit======
    Qneg[2] = c[2][0]*Q[0] + c[2][1]*Q[1] + c[2][2]*Q[2];
    Qneg[1] = c[1][0]*Q[1] + c[1][1]*Q[2] + c[1][2]*Q[3];
    Qneg[0] = c[0][0]*Q[2] + c[0][1]*Q[3] + c[0][2]*Q[4];

    beta[0] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[1] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    beta[2] = (Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[2] - 4.0*Q[3]     + Q[4])*(3.0*Q[2] - 4.0*Q[3]     + Q[4]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[1]                - Q[3])*(    Q[1]                - Q[3]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[0] - 4.0*Q[1] + 3.0*Q[2])*(    Q[0] - 4.0*Q[1] + 3.0*Q[2]);
    alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dneg[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[6] = Qneg[0]*omg[0] + Qneg[1]*omg[1] + Qneg[2]*omg[2];
}



void local_WENO_5_interright(double h, double Q[8])
{
  double d[3][3], eps = 1e-6;
  d[1][0] = 1.0/3.0;
  d[1][1] = 5.0/6.0;
  d[1][2] = -1.0/6.0;
  d[2][0] = d[1][2];
  d[2][1] = d[1][1];
  d[2][2] = d[1][0];
  d[0][0] = 11.0/6.0;
  d[0][1] = -7.0/6.0;
  d[0][2] = 1.0/3.0;

  double Dpos[3] = {0.1, 0.6, 0.3};//dL

  double alp[3], omg[3], beta[3], sum;
  double Qpos[3], dQpos[3];


  //=====Reconstruction for The Positive Limit======
    Qpos[2] = d[2][0]*Q[1] + d[2][1]*Q[2] + d[2][2]*Q[3];
    Qpos[1] = d[1][0]*Q[2] + d[1][1]*Q[3] + d[1][2]*Q[4];
    Qpos[0] = d[0][0]*Q[3] + d[0][1]*Q[4] + d[0][2]*Q[5];

    beta[0] = (Q[3] - 2.0*Q[4] + Q[5])*(Q[3] - 2.0*Q[4] + Q[5]);
    beta[1] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[2] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[3] - 4.0*Q[4]     + Q[5])*(3.0*Q[3] - 4.0*Q[4]     + Q[5]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[2]                - Q[4])*(    Q[2]                - Q[4]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[1] - 4.0*Q[2] + 3.0*Q[3])*(    Q[1] - 4.0*Q[2] + 3.0*Q[3]);
    alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dpos[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[7] = Qpos[0]*omg[0] + Qpos[1]*omg[1] + Qpos[2]*omg[2];
}



void local_WENO_5_inter(double h, double Q[8])
{
  double c[3][3], d[3][3], e[3][3], f[3][3], eps = 1e-6;
  c[0][0] = 1.0/3.0;
  c[0][1] = 5.0/6.0;
  c[0][2] = -1.0/6.0;
  c[1][0] = c[0][2];
  c[1][1] = c[0][1];
  c[1][2] = c[0][0];
  c[2][0] = 1.0/3.0;
  c[2][1] = -7.0/6.0;
  c[2][2] = 11.0/6.0;
  d[1][0] = c[0][0];
  d[1][1] = c[0][1];
  d[1][2] = c[0][2];
  d[2][0] = c[1][0];
  d[2][1] = c[1][1];
  d[2][2] = c[1][2];
  d[0][0] = 11.0/6.0;
  d[0][1] = -7.0/6.0;
  d[0][2] = 1.0/3.0;

  double Dneg[3] = {0.3, 0.6, 0.1};//dR
  double Dpos[3] = {0.1, 0.6, 0.3};//dL

  double alp[3], omg[3], beta[3], sum;
  double Qneg[3], Qpos[3];


  //=====Reconstruction for The Negative Limit======
    Qneg[2] = c[2][0]*Q[0] + c[2][1]*Q[1] + c[2][2]*Q[2];
    Qneg[1] = c[1][0]*Q[1] + c[1][1]*Q[2] + c[1][2]*Q[3];
    Qneg[0] = c[0][0]*Q[2] + c[0][1]*Q[3] + c[0][2]*Q[4];

    beta[0] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[1] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    beta[2] = (Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[2] - 4.0*Q[3]     + Q[4])*(3.0*Q[2] - 4.0*Q[3]     + Q[4]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[1]                - Q[3])*(    Q[1]                - Q[3]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[0] - 4.0*Q[1] + 3.0*Q[2])*(    Q[0] - 4.0*Q[1] + 3.0*Q[2]);
    alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dneg[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[6] = Qneg[0]*omg[0] + Qneg[1]*omg[1] + Qneg[2]*omg[2];


  //=====Reconstruction for The Positive Limit======
    Qpos[2] = d[2][0]*Q[1] + d[2][1]*Q[2] + d[2][2]*Q[3];
    Qpos[1] = d[1][0]*Q[2] + d[1][1]*Q[3] + d[1][2]*Q[4];
    Qpos[0] = d[0][0]*Q[3] + d[0][1]*Q[4] + d[0][2]*Q[5];

    beta[0] = (Q[3] - 2.0*Q[4] + Q[5])*(Q[3] - 2.0*Q[4] + Q[5]);
    beta[1] = (Q[2] - 2.0*Q[3] + Q[4])*(Q[2] - 2.0*Q[3] + Q[4]);
    beta[2] = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
    beta[0] = (13.0/12.0)*beta[0] + 0.25*(3.0*Q[3] - 4.0*Q[4]     + Q[5])*(3.0*Q[3] - 4.0*Q[4]     + Q[5]);
    beta[1] = (13.0/12.0)*beta[1] + 0.25*(    Q[2]                - Q[4])*(    Q[2]                - Q[4]);
    beta[2] = (13.0/12.0)*beta[2] + 0.25*(    Q[1] - 4.0*Q[2] + 3.0*Q[3])*(    Q[1] - 4.0*Q[2] + 3.0*Q[3]);
    alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
    alp[2] = Dpos[2] / ((eps+beta[2])*(eps+beta[2]));
    sum = alp[0] + alp[1] + alp[2];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    omg[2] = alp[2] / sum;
    Q[7] = Qpos[0]*omg[0] + Qpos[1]*omg[1] + Qpos[2]*omg[2];
}
