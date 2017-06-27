#include <math.h>


void flux_RF(double const running_info[], int const m, double const h, double const gamma,
	     double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[])
{
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];
  int j, k, a;

  double rho_L, rho_R, u_L, u_R, p_L, p_R, c_L, c_R;
  double rho_star_L, rho_star_R, u_star_L, u_star_R, p_star_L, p_star_R, c_star_L, c_star_R;
  double alp1, alp3, add, beta1, beta3;

  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double f1[m+6], f2[m+6], f3[m+6];
  double Q1[8], Q2[8], Q3[8];
  double u, p, g1, g2, g3;

  double H_star, u_star, c_star;
  double H[m+6];
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
    u = W2[j]/W1[j];
    p = (W3[j] - 0.5*W2[j]*u)*(gamma-1.0);
    f1[j] = W2[j];
    f2[j] = W2[j]*u + p;
    f3[j] = (W3[j]+p)*u;
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
      Q1[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((f2[j-3+k]-0.5*f1[j-3+k]*u_star)*u_star - f3[j-3+k]) / c_star/c_star + f1[j-3+k];
    }


  //=====Reconstruction & Recomposition========
    rho_L = W1[j-1];
      u_L = W2[j-1]/rho_L;
      p_L = (gamma-1.0)*(W3[j-1] - 0.5*W2[j-1]*u_L);
      c_L = sqrt(gamma*p_L/rho_L);
    if(u_L-c_L < 0.0)//possible left transonic rarefaction
    {
      alp1 = (Q1[3]-Q1[2])/(u_star-c_star);
      /*
       *   Q1[3] - Q1[2] = beta^1_R - beta^1_L
       * = dF.L1
       * = (A dU).L1
       * = dU^t A^t L1
       * = dU^t [L1 L2 L3] Lambda R L1
       * = [alp1 alp2 alp3] Lambda [1 0 0]^t
       * = [lamb1*alp1 lamb2*alp2 lamb3*alp3] [1 0 0]^t
       * = lamb1*alp1
       * ==> alp1 = (Q1[3]-Q1[2])/lamb1
       */
      rho_star_L = W1[j-1] + alp1;
      u_star_L = (W2[j-1] + alp1*(u_star-c_star));
      u_star_L = u_star_L / rho_star_L;
      p_star_L = W3[j-1] + alp1*(H_star-u_star*c_star);
      p_star_L = (gamma-1.0)*(p_star_L - 0.5*rho_star_L*u_star_L*u_star_L);
      c_star_L = sqrt(gamma*p_star_L/rho_star_L);
      /*  ^    U*L = UL + alp1 * R1
       *  |  compute u_star_L-c_star_L to decide whether there is a left
       *     transonic rarefaction
       */
      if(u_star_L-c_star_L > 0.0)//left transonic rarefaction
      {
	local_WENO_5_interleft_Z(h, Q1);
	local_WENO_5_interleft_Z(h, Q2);
	local_WENO_5_interleft_Z(h, Q3);
	g1 = Q1[6];
	g2 = Q2[6];
	g3 = Q3[6];
	g1+= alp1*(u_L-c_L)*((u_star_L-c_star_L) - (u_star-c_star))/((u_star_L-c_star_L) - (u_L-c_L));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }
    rho_R = W1[j];
      u_R = W2[j]/rho_R;
      p_R = (gamma-1.0)*(W3[j] - 0.5*W2[j]*u_R);
      c_R = sqrt(gamma*p_R/rho_R);
    if(u_R+c_R > 0.0)//possible right transonic rarefaction
    {
      alp3 = (Q3[3]-Q3[2])/(u_star+c_star);
      rho_star_R = W1[j] - alp3;
      u_star_R = (W2[j] - alp3*(u_star+c_star));
      u_star_R = u_star_R / rho_star_R;
      p_star_R = W3[j] - alp3*(H_star+u_star*c_star);
      p_star_R = (gamma-1.0)*(p_star_R - 0.5*rho_star_R*u_star_R*u_star_R);
      c_star_R = sqrt(gamma*p_star_R/rho_star_R);
      /*  ^    U*R = UR - alp3 * R3
       *  |  compute u_star_R+c_star_R to decide whether there is a right
       *     transonic rarefaction
       */
      if(u_star_R+c_star_R < 0.0)//right transonic rarefaction
      {
	local_WENO_5_interright_Z(h, Q1);
	local_WENO_5_interright_Z(h, Q2);
	local_WENO_5_interright_Z(h, Q3);
	g1 = Q1[7];
	g2 = Q2[7];
	g3 = Q3[7];
	g3-= alp3*(u_R+c_R)*((u_star+c_star) - (u_star_R+c_star_R))/((u_R+c_R) - (u_star_R+c_star_R));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }

    //non transonic rarefaction
    if(u_star-c_star > 0.0)
    {
      local_WENO_5_interleft_Z(h, Q1);
      local_WENO_5_interleft_Z(h, Q2);
      local_WENO_5_interleft_Z(h, Q3);
      g1 = Q1[6];
      g2 = Q2[6];
      g3 = Q3[6];
    }
    else if(u_star > 0.0)
    {
      local_WENO_5_interright_Z(h, Q1);
      local_WENO_5_interleft_Z(h, Q2);
      local_WENO_5_interleft_Z(h, Q3);
      g1 = Q1[7];
      g2 = Q2[6];
      g3 = Q3[6];
    }
    else if(u_star+c_star > 0.0)
    {
      local_WENO_5_interright_Z(h, Q1);
      local_WENO_5_interright_Z(h, Q2);
      local_WENO_5_interleft_Z(h, Q3);
      g1 = Q1[7];
      g2 = Q2[7];
      g3 = Q3[6];
    }
    else
    {
      local_WENO_5_interright_Z(h, Q1);
      local_WENO_5_interright_Z(h, Q2);
      local_WENO_5_interright_Z(h, Q3);
      g1 = Q1[7];
      g2 = Q2[7];
      g3 = Q3[7];
    }

    F1[j-3] = g1+g2+g3;
    F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
    F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
    a = 0;
  }

  if(bod < 0)
  {
    F1[0] = 0.0;
    F1[m] = 0.0;
    F3[0] = 0.0;
    F3[m] = 0.0;
  }
}


/*
 * We use this dual Roe solver to solve dual Riemann problem
 *     u_t - f(u)_x = 0
 * Its Roe linearization is
 *     u_t - A u_x = 0
 * The eigenvector of -A is the same as those of A. The eigenvalues
 * of -A is the opposite signs, i.e. they are
 *     [eta_3, eta_2, eta_1] = [-u-c, -u, -u+c]
 * corresponding to the right eigenvectors
 *                [R3, R2, R1].
 * More details can be seen in our document.
 */
void flux_RF_dual(double const running_info[], int const m, double const h, double const gamma,
		  double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[])
{
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];
  int j, k, a;

  double rho_L, rho_R, u_L, u_R, p_L, p_R, c_L, c_R;
  double rho_star_L, rho_star_R, u_star_L, u_star_R, p_star_L, p_star_R, c_star_L, c_star_R;
  double alp1, alp3, add, beta1, beta3;

  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double f1[m+6], f2[m+6], f3[m+6];
  double Q1[8], Q2[8], Q3[8];
  double u, p, g1, g2, g3;

  double H_star, u_star, c_star;
  double H[m+6];
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
    u = W2[j]/W1[j];
    p = (W3[j] - 0.5*W2[j]*u)*(gamma-1.0);
    f1[j] = W2[j];
    f2[j] = W2[j]*u + p;
    f3[j] = (W3[j]+p)*u;
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
      Q1[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((f2[j-3+k]-0.5*f1[j-3+k]*u_star)*u_star - f3[j-3+k]) / c_star/c_star + f1[j-3+k];
    }


  //=====Reconstruction & Recomposition========
    rho_L = W1[j-1];
    u_L = W2[j-1]/rho_L;
    p_L = (gamma-1.0)*(W3[j-1] - 0.5*W2[j-1]*u_L);
    c_L = sqrt(gamma*p_L/rho_L);
    if(-u_L-c_L < 0.0)//possible left transonic rarefaction
    {
      alp3 = (Q3[3]-Q3[2])/(u_star+c_star);
      rho_star_L = W1[j] + alp3;
      u_star_L = (W2[j] + alp3*(u_star+c_star));
      p_star_L = W3[j] + alp3*(H_star+u_star*c_star);
      u_star_L = u_star_L / rho_star_L;
      p_star_L = (gamma-1.0)*(p_star_L - 0.5*rho_star_L*u_star_L*u_star_L);
      c_star_L = sqrt(gamma*p_star_L/rho_star_L);
      /*  ^    U*L = UL + alp3 * R3
       *  |  compute -u_star_L-c_star_L to decide whether there is a left
       *     transonic rarefaction
       */
      if(-u_star_L-c_star_L > 0.0)//left transonic rarefaction
      {
	local_WENO_5_interleft_Z(h, Q1);
	local_WENO_5_interleft_Z(h, Q2);
	local_WENO_5_interleft_Z(h, Q3);
	g1 = Q1[6];
	g2 = Q2[6];
	g3 = Q3[6];
	g3+= alp3*(-u_L-c_L)*((-u_star_L-c_star_L) - (-u_star-c_star))/((-u_star_L-c_star_L) - (-u_L-c_L));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }
    rho_R = W1[j];
    u_R = W2[j]/rho_R;
    p_R = (gamma-1.0)*(W3[j] - 0.5*W2[j]*u_R);
    c_R = sqrt(gamma*p_R/rho_R);
    if(-u_R+c_R > 0.0)//possible right transonic rarefaction
    {
      alp1 = (Q1[3]-Q1[2])/(u_star-c_star);
      /*
       *   Q1[3] - Q1[2] = beta^1_R - beta^1_L
       * = dF.L1
       * = (A dU).L1
       * = dU^t A^t L1
       * = dU^t [L1 L2 L3] Lambda R L1
       * = [alp1 alp2 alp3] Lambda [1 0 0]^t
       * = [lamb1*alp1 lamb2*alp2 lamb3*alp3] [1 0 0]^t
       * = lamb1*alp1
       * ==> alp1 = (Q1[3]-Q1[2])/lamb1
       */
      rho_star_R = W1[j-1] - alp1;
      u_star_R = (W2[j-1] - alp1*(u_star-c_star));
      p_star_R = W3[j-1] - alp1*(H_star-u_star*c_star);
      u_star_R = u_star_R / rho_star_R;
      p_star_R = (gamma-1.0)*(p_star_R - 0.5*rho_star_R*u_star_R*u_star_R);
      c_star_R = sqrt(gamma*p_star_R/rho_star_R);
      /*  ^    U*R = UR - alp1 * R1
       *  |  compute u_star_R-c_star_R to decide whether there is a right
       *     transonic rarefaction
       */
      if(-u_star_R+c_star_R < 0.0)//right transonic rarefaction
      {
	local_WENO_5_interright_Z(h, Q1);
	local_WENO_5_interright_Z(h, Q2);
	local_WENO_5_interright_Z(h, Q3);
	g1 = Q1[7];
	g2 = Q2[7];
	g3 = Q3[7];
	g1-= alp1*(-u_R+c_R)*((-u_star+c_star) - (-u_star_R+c_star_R))/((-u_R+c_R) - (-u_star_R+c_star_R));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }

    //non transonic rarefaction
    if(u_star-c_star > 0.0)
    {
      local_WENO_5_interright_Z(h, Q1);
      local_WENO_5_interright_Z(h, Q2);
      local_WENO_5_interright_Z(h, Q3);
      g1 = Q1[7];
      g2 = Q2[7];
      g3 = Q3[7];
    }
    else if(u_star > 0.0)
    {
      local_WENO_5_interleft_Z(h, Q1);
      local_WENO_5_interright_Z(h, Q2);
      local_WENO_5_interright_Z(h, Q3);
      g1 = Q1[6];
      g2 = Q2[7];
      g3 = Q3[7];
    }
    else if(u_star+c_star > 0.0)
    {
      local_WENO_5_interleft_Z(h, Q1);
      local_WENO_5_interleft_Z(h, Q2);
      local_WENO_5_interright_Z(h, Q3);
      g1 = Q1[6];
      g2 = Q2[6];
      g3 = Q3[7];
    }
    else
    {
      local_WENO_5_interleft_Z(h, Q1);
      local_WENO_5_interleft_Z(h, Q2);
      local_WENO_5_interleft_Z(h, Q3);
      g1 = Q1[6];
      g2 = Q2[6];
      g3 = Q3[6];
    }

    F1[j-3] = g1+g2+g3;
    F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
    F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
    a = 0;
  }

  if(bod < 0)
  {
    F1[0] = 0.0;
    F1[m] = 0.0;
    F3[0] = 0.0;
    F3[m] = 0.0;
  }
}



void flux_RF_1st(int const running_info[], int const m, double const h, double const gamma,
	     double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[])
{
  int j, k, a;
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];


  double rho_L, rho_R, u_L, u_R, p_L, p_R, c_L, c_R;
  double rho_star_L, rho_star_R, u_star_L, u_star_R, p_star_L, p_star_R, c_star_L, c_star_R;
  double alp1, alp3, add, beta1, beta3;
  double QQ[3][2];

  double W1[m+6], W2[m+6], W3[m+6]; //W1=rho, W2=mom, W3=ene
  double f1[m+6], f2[m+6], f3[m+6];
  double Q1[6], Q2[6], Q3[6];
  double u, p, g1, g2, g3;

  double H_star, u_star, c_star;
  double H[m+6];
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
    u = W2[j]/W1[j];
    p = (W3[j] - 0.5*W2[j]*u)*(gamma-1.0);
    f1[j] = W2[j];
    f2[j] = W2[j]*u + p;
    f3[j] = (W3[j]+p)*u;
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
      Q1[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star + 1.0);
      Q1[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star - 1.0);
      Q1[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q1[k] = Q1[k] / c_star;
      Q3[k] = f1[j-3+k] * 0.5*u_star* (0.5*(gamma-1.0)*u_star/c_star - 1.0);
      Q3[k]+= f2[j-3+k] * 0.5 * ((1.0-gamma)*u_star/c_star + 1.0);
      Q3[k]+= f3[j-3+k] * 0.5 * (gamma-1.0) / c_star;
      Q3[k] = Q3[k] / c_star;
      Q2[k] = (gamma-1.0) * ((f2[j-3+k]-0.5*f1[j-3+k]*u_star)*u_star - f3[j-3+k]) / c_star/c_star + f1[j-3+k];
    }


  //=====Reconstruction & Recomposition========
    rho_L = W1[j-1];
      u_L = W2[j-1]/rho_L;
      p_L = (gamma-1.0)*(W3[j-1] - 0.5*W2[j-1]*u_L);
      c_L = sqrt(gamma*p_L/rho_L);
    if(u_L-c_L < 0.0)//possible left transonic rarefaction
    {
      alp1 = (Q1[3]-Q1[2])/(u_star-c_star);
      rho_star_L = W1[j-1] + alp1;
        u_star_L = (W2[j-1] + alp1*(u_star-c_star)) / rho_star_L;
        p_star_L = W3[j-1] + alp1*(H_star-u_star*c_star);
        p_star_L = (gamma-1.0)*(p_star_L - 0.5*rho_star_L*u_star_L*u_star_L);
	c_star_L = sqrt(gamma*p_star_L/rho_star_L);
      if(u_star_L-c_star_L > 0.0)//left transonic rarefaction
      {
	//local_WENO_5_d_interleft(h, Q1, QQ[0]);
	//local_WENO_5_d_interleft(h, Q2, QQ[1]);
	//local_WENO_5_d_interleft(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
	g1 = QQ[0][0];
	g2 = QQ[1][0];
	g3 = QQ[2][0];
	g1+= alp1*(u_L-c_L)*((u_star_L-c_star_L) - (u_star-c_star))/((u_star_L-c_star_L) - (u_L-c_L));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }
    rho_R = W1[j];
      u_R = W2[j]/rho_R;
      p_R = (gamma-1.0)*(W3[j] - 0.5*W2[j]*u_R);
      c_R = sqrt(gamma*p_R/rho_R);
    if(u_R+c_R > 0.0)//possible right transonic rarefaction
    {
      alp3 = (Q3[3]-Q3[2])/(u_star+c_star);
      rho_star_R = W1[j] - alp3;
        u_star_R = (W2[j] - alp3*(u_star+c_star)) / rho_star_R;
        p_star_R = W3[j] - alp3*(H_star+u_star*c_star);
        p_star_R = (gamma-1.0)*(p_star_R - 0.5*rho_star_R*u_star_R*u_star_R);
	c_star_R = sqrt(gamma*p_star_R/rho_star_R);
      if(u_star_R+c_star_R < 0.0)//right transonic rarefaction
      {
	//local_WENO_5_d_interright(h, Q1, QQ[0]);
	//local_WENO_5_d_interright(h, Q2, QQ[1]);
	//local_WENO_5_d_interright(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
	g1 = QQ[0][1];
	g2 = QQ[1][1];
	g3 = QQ[2][1];
	g3-= alp3*(u_R+c_R)*((u_star+c_star) - (u_star_R+c_star_R))/((u_R+c_R) - (u_star_R+c_star_R));
	F1[j-3] = g1+g2+g3;
	F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
	F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
	continue;
      }
    }

    //non transonic rarefaction
    if(u_star-c_star > 0.0)
    {
      //local_WENO_5_d_interleft(h, Q1, QQ[0]);
      //local_WENO_5_d_interleft(h, Q2, QQ[1]);
      //local_WENO_5_d_interleft(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
      g1 = QQ[0][0];
      g2 = QQ[1][0];
      g3 = QQ[2][0];
    }
    else if(u_star > 0.0)
    {
      //local_WENO_5_d_interright(h, Q1, QQ[0]);
      //local_WENO_5_d_interleft(h, Q2, QQ[1]);
      //local_WENO_5_d_interleft(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
      g1 = QQ[0][1];
      g2 = QQ[1][0];
      g3 = QQ[2][0];
    }
    else if(u_star+c_star > 0.0)
    {
      //local_WENO_5_d_interright(h, Q1, QQ[0]);
      //local_WENO_5_d_interright(h, Q2, QQ[1]);
      //local_WENO_5_d_interleft(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
      g1 = QQ[0][1];
      g2 = QQ[1][1];
      g3 = QQ[2][0];
    }
    else
    {
      //local_WENO_5_d_interright(h, Q1, QQ[0]);
      //local_WENO_5_d_interright(h, Q2, QQ[1]);
      //local_WENO_5_d_interright(h, Q3, QQ[2]);
	QQ[0][0] = Q1[2];  QQ[0][1] = Q1[3];
	QQ[1][0] = Q2[2];  QQ[1][1] = Q2[3];
	QQ[2][0] = Q3[2];  QQ[2][1] = Q3[3];
      g1 = QQ[0][1];
      g2 = QQ[1][1];
      g3 = QQ[2][1];
    }

    F1[j-3] = g1+g2+g3;
    F2[j-3] = u_star*F1[j-3] + c_star*(g3 - g1);
    F3[j-3] = H_star*(g1+g3) + u_star*c_star*(g3-g1) + 0.5*u_star*u_star*g2;
    a = 0;
  }

  if(bod < 0)
  {
    F1[0] = 0.0;
    F1[m] = 0.0;
    F3[0] = 0.0;
    F3[m] = 0.0;
  }
}
