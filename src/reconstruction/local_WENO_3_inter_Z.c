#include <math.h>


void local_WENO_3_inter_Z(double h, double Q[6])
{
  double c[2][2], d[2][2], eps = 1e-40;
  eps = 1e-6;
  c[0][0] = 0.5;
  c[0][1] = 0.5;
  c[1][0] = -0.5;
  c[1][1] = 1.5;
  d[0][0] = 1.5;
  d[0][1] = -0.5;
  d[1][0] = 0.5;
  d[1][1] = 0.5;

  double Dneg[2] = {2.0/3.0, 1.0/3.0};//dR
  double Dpos[2] = {1.0/3.0, 2.0/3.0};//dL

  double alp[2], omg[2], beta[2], tau, sum;
  double Qneg[2], Qpos[2]; 

  /*
   * Stencil
   *   #==0==#==1==#==2==#
   *   |----[1]----|
   *         |----[0]----|
   */

  //=====Reconstruction for The Negative Limit======
    Qneg[1] = c[1][0]*Q[0] + c[1][1]*Q[1];
    Qneg[0] = c[0][0]*Q[1] + c[0][1]*Q[2];

    beta[0] = (Q[2] - Q[1])*(Q[2] - Q[1]);
    beta[1] = (Q[1] - Q[0])*(Q[1] - Q[0]);
    /*
    tau = fabs(beta[1] - beta[0]);
    beta[0] = tau / (beta[0] + eps);
    beta[1] = tau / (beta[1] + eps);

    alp[0] = Dneg[0] * (1 + beta[0]*beta[0]);
    alp[1] = Dneg[1] * (1 + beta[1]*beta[1]);
    /*/
    alp[0] = Dneg[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dneg[1] / ((eps+beta[1])*(eps+beta[1]));
    //*/
    sum = alp[0] + alp[1];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    Q[4] = Qneg[0]*omg[0] + Qneg[1]*omg[1];
    //Q[4] = Qneg[0]*Dneg[0] + Qneg[1]*Dneg[1];
    //Q[4] = (-Q[0] + 5.0*Q[1] + 2.0*Q[2])/6.0;


  /*
   * Stencil
   *   #==1==#==2==#==3==#
   *   |----[1]----|
   *         |----[0]----|
   */

  //=====Reconstruction for The Positive Limit======
    Qpos[1] = d[1][0]*Q[1] + d[1][1]*Q[2];
    Qpos[0] = d[0][0]*Q[2] + d[0][1]*Q[3];

    beta[0] = (Q[2] - Q[3])*(Q[2] - Q[3]);
    beta[1] = (Q[2] - Q[1])*(Q[2] - Q[1]);
    /*
    tau = fabs(beta[1] - beta[0]);
    beta[0] = tau / (beta[0] + eps);
    beta[1] = tau / (beta[1] + eps);

    alp[0] = Dpos[0] * (1 + beta[0]*beta[0]);
    alp[1] = Dpos[1] * (1 + beta[1]*beta[1]);
    /*/
    alp[0] = Dpos[0] / ((eps+beta[0])*(eps+beta[0]));
    alp[1] = Dpos[1] / ((eps+beta[1])*(eps+beta[1]));
    //*/
    sum = alp[0] + alp[1];
    omg[0] = alp[0] / sum;
    omg[1] = alp[1] / sum;
    Q[5] = Qpos[0]*omg[0] + Qpos[1]*omg[1];
    //Q[5] = Qpos[0]*Dpos[0] + Qpos[1]*Dpos[1];
    //Q[5] = (2.0*Q[1] + 5.0*Q[2] - Q[3])/6.0;
}
