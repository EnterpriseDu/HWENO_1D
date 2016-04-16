#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void local_HWENO_5_inter_d(double h, double Q[6], double DQ[6])
{
  double eps = 1e-6;
  double Dneg[3] = {0.1125, 0.525, 0.3625};
  double Dpos[3] = {0.525, 0.1125, 0.3625};
  double dDneg[3] = {1.0/18.0, 5.0/6.0, 1.0/9.0};
  double dDpos[3] = {5.0/6.0, 1.0/18.0, 1.0/9.0};

  double alp[3], omg[3], beta[3], sum;
  double Qneg[3], Qpos[3], dQneg[3], dQpos[3];


  /*
   * Stencil
   *   #=====#=====#=====#
   *   |----[0]----|
   *         |----[1]----|
   *   |-------[2]-------|
   */

  //=====Reconstruction for The Negative Limit======
  Qneg[0] = (-7.0*Q[0] + 13.0*Q[1] - 4.0*DQ[0]) / 6.0;
  Qneg[1] = (Q[1] + 5.0*Q[2] - 2.0*DQ[2]) / 6.0;
  Qneg[2] = (-Q[0] + 5.0*Q[1] + 2.0*Q[2]) / 6.0;

  dQneg[0] =  (4.0* (Q[0] - Q[1]) + 1.5*DQ[0] + 3.5*DQ[1] )/h;
  dQneg[1] =  (2.0* (Q[2] - Q[1]) - 0.5*DQ[1]  - 0.5*DQ[2])/h;
  dQneg[2] =  (0.25*(Q[0] - 4.0*Q[1] + 3.0*Q[2]) + 0.5*DQ[1])/h;

  beta[0]  = 4.0*(3.0*(Q[0] - Q[1]) + DQ[0] + 2.0*DQ[1])*(3.0*(Q[0] - Q[1]) + DQ[0] + 2.0*DQ[1]);
  beta[1]  = 4.0*(3.0*(Q[2] - Q[1]) - DQ[2] - 2.0*DQ[1])*(3.0*(Q[2] - Q[1]) - DQ[2] - 2.0*DQ[1]);
  beta[0] += 9.75*( 2.0*(Q[0] - Q[1]) + DQ[0] + DQ[1])*( 2.0*(Q[0] - Q[1]) + DQ[0] + DQ[1]);
  beta[1] += 9.75*(-2.0*(Q[2] - Q[1]) + DQ[2] + DQ[1])*(-2.0*(Q[2] - Q[1]) + DQ[2] + DQ[1]);
  beta[2]  = (Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
  beta[2] += 2.4375*(Q[2] - Q[0] - 2.0*DQ[1])*(Q[2] - Q[0] - 2.0*DQ[1]);
  alp[0] = dDneg[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = dDneg[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = dDneg[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  DQ[4] = dQneg[0]*omg[0] + dQneg[1]*omg[1] + dQneg[2]*omg[2];

  beta[0] =        (-2.0*Q[0] + 2.0*Q[1] - DQ[0])*(-2.0*Q[0] + 2.0*Q[1] - DQ[0]);
  beta[0]+= (13.0/3.0)*(-Q[0] +     Q[1] - DQ[0])*(    -Q[0] +     Q[1] - DQ[0]);
  beta[1] =         (2.0*Q[2] - 2.0*Q[1] - DQ[2])*( 2.0*Q[2] - 2.0*Q[1] - DQ[2]);
  beta[1]+= (13.0/3.0)*( Q[2] -     Q[1] - DQ[2])*(     Q[2] -     Q[1] - DQ[2]);
  beta[2] = (13.0/12.0)*(Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
  beta[2]+= 0.25*(Q[0] - Q[2])*(Q[0] - Q[2]);
  alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = Dneg[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  Q[4] = Qneg[0]*omg[0] + Qneg[1]*omg[1] + Qneg[2]*omg[2];



  //=====Reconstruction for The Positive Limit======
  Qpos[1] = (-7.0*Q[3] + 13.0*Q[2] + 4.0*DQ[3]) / 6.0;
  Qpos[0] = (Q[2] + 5.0*Q[1] + 2.0*DQ[1]) / 6.0;
  Qpos[2] = (-Q[3] + 5.0*Q[2] + 2.0*Q[1]) / 6.0;

  dQpos[1] = -(4.0* (Q[3] - Q[2]) - 1.5*DQ[3] - 3.5*DQ[2] )/h;
  dQpos[0] = -(2.0* (Q[1] - Q[2]) + 0.5*DQ[2]  + 0.5*DQ[1])/h;
  dQpos[2] = -(0.25*(Q[3] - 4.0*Q[2] + 3.0*Q[1]) - 0.5*DQ[2])/h;

  beta[0]  = 4.0*(3.0*(Q[1] - Q[2]) + DQ[1] + 2.0*DQ[2])*(3.0*(Q[1] - Q[2]) + DQ[1] + 2.0*DQ[2]);
  beta[1]  = 4.0*(3.0*(Q[3] - Q[2]) - DQ[3] - 2.0*DQ[2])*(3.0*(Q[3] - Q[2]) - DQ[3] - 2.0*DQ[2]);
  beta[0] += 9.75*( 2.0*(Q[1] - Q[2]) + DQ[1] + DQ[2])*( 2.0*(Q[1] - Q[2]) + DQ[1] + DQ[2]);
  beta[1] += 9.75*(-2.0*(Q[3] - Q[2]) + DQ[3] + DQ[2])*(-2.0*(Q[3] - Q[2]) + DQ[3] + DQ[2]);
  beta[2]  = (Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
  beta[2] += 2.4375*(Q[3] - Q[1] - 2.0*DQ[2])*(Q[3] - Q[1] - 2.0*DQ[2]);
  alp[0] = dDpos[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = dDpos[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = dDpos[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  DQ[5] = dQpos[0]*omg[0] + dQpos[1]*omg[1] + dQpos[2]*omg[2];

  beta[0] =        (-2.0*Q[1] + 2.0*Q[2] - DQ[1])*(-2.0*Q[1] + 2.0*Q[2] - DQ[1]);
  beta[0]+= (13.0/3.0)*(-Q[1] +     Q[2] - DQ[1])*(    -Q[1] +     Q[2] - DQ[1]);
  beta[1] =         (2.0*Q[3] - 2.0*Q[2] - DQ[3])*( 2.0*Q[3] - 2.0*Q[2] - DQ[3]);
  beta[1]+= (13.0/3.0)*( Q[3] -     Q[2] - DQ[3])*(     Q[3] -     Q[2] - DQ[3]);
  beta[2] = (13.0/12.0)*(Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
  beta[2]+= 0.25*(Q[1] - Q[3])*(Q[1] - Q[3]);
  alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = Dpos[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  Q[5] = Qpos[0]*omg[0] + Qpos[1]*omg[1] + Qpos[2]*omg[2];
}



void local_HWENO_5_inter(double h, double Q[6], double DQ[4])
{
  double eps = 1e-40;
  double Dneg[3] = {0.1125, 0.525, 0.3625};
  double Dpos[3] = {0.525, 0.1125, 0.3625};

  double alp[3], omg[3], beta[3], sum;
  double Qneg[3], Qpos[3];


  /*
   * Stencil
   *   #=====#=====#=====#
   *   |----[0]----|
   *         |----[1]----|
   *   |-------[2]-------|
   */

  //=====Reconstruction for The Negative Limit======
  Qneg[0] = (-7.0*Q[0] + 13.0*Q[1] - 4.0*DQ[0]) / 6.0;
  Qneg[1] = (Q[1] + 5.0*Q[2] - 2.0*DQ[2]) / 6.0;
  Qneg[2] = (-Q[0] + 5.0*Q[1] + 2.0*Q[2]) / 6.0;

  beta[0] =        (-2.0*Q[0] + 2.0*Q[1] - DQ[0])*(-2.0*Q[0] + 2.0*Q[1] - DQ[0]);
  beta[0]+= (13.0/3.0)*(-Q[0] +     Q[1] - DQ[0])*(    -Q[0] +     Q[1] - DQ[0]);
  beta[1] =         (2.0*Q[2] - 2.0*Q[1] - DQ[2])*( 2.0*Q[2] - 2.0*Q[1] - DQ[2]);
  beta[1]+= (13.0/3.0)*( Q[2] -     Q[1] - DQ[2])*(     Q[2] -     Q[1] - DQ[2]);
  beta[2] = (13.0/12.0)*(Q[0] - 2.0*Q[1] + Q[2])*(Q[0] - 2.0*Q[1] + Q[2]);
  beta[2]+= 0.25*(Q[0] - Q[2])*(Q[0] - Q[2]);
  alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = Dneg[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  Q[4] = Qneg[0]*omg[0] + Qneg[1]*omg[1] + Qneg[2]*omg[2];



  //=====Reconstruction for The Positive Limit======
  Qpos[1] = (-7.0*Q[3] + 13.0*Q[2] + 4.0*DQ[3]) / 6.0;
  Qpos[0] = (Q[2] + 5.0*Q[1] + 2.0*DQ[1]) / 6.0;
  Qpos[2] = (-Q[3] + 5.0*Q[2] + 2.0*Q[1]) / 6.0;

  beta[0] =        (-2.0*Q[1] + 2.0*Q[2] - DQ[1])*(-2.0*Q[1] + 2.0*Q[2] - DQ[1]);
  beta[0]+= (13.0/3.0)*(-Q[1] +     Q[2] - DQ[1])*(    -Q[1] +     Q[2] - DQ[1]);
  beta[1] =         (2.0*Q[3] - 2.0*Q[2] - DQ[3])*( 2.0*Q[3] - 2.0*Q[2] - DQ[3]);
  beta[1]+= (13.0/3.0)*( Q[3] -     Q[2] - DQ[3])*(     Q[3] -     Q[2] - DQ[3]);
  beta[2] = (13.0/12.0)*(Q[1] - 2.0*Q[2] + Q[3])*(Q[1] - 2.0*Q[2] + Q[3]);
  beta[2]+= 0.25*(Q[1] - Q[3])*(Q[1] - Q[3]);
  alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
  alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
  alp[2] = Dpos[2] / ((eps+beta[2])*(eps+beta[2]));
  sum = alp[0] + alp[1] + alp[2];
  omg[0] = alp[0] / sum;
  omg[1] = alp[1] / sum;
  omg[2] = alp[2] / sum;
  Q[5] = Qpos[0]*omg[0] + Qpos[1]*omg[1] + Qpos[2]*omg[2];
}
