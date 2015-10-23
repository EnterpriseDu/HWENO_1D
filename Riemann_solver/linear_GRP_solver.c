#include <math.h>
#include <stdio.h>





void linear_GRP_solver
(double wave_speed[2], double D[4], double U[4], double lambda, double gamma, double eps,
 double rho_L, double u_L, double v_L, double p_L,
 double rho_R, double u_R, double v_R, double p_R,
 double d_rho_L, double d_u_L, double d_v_L, double d_p_L,
 double d_rho_R, double d_u_R, double d_v_R, double d_p_R,
 double t_rho_L, double t_u_L, double t_v_L, double t_p_L,
 double t_rho_R, double t_u_R, double t_v_R, double t_p_R)
{
  double dist;
  double c_L, c_R, c_frac, cricit, c_star, c_square;
  int CRW[2];
  double d_Phi, d_Psi, TdS, VAR, RIEL, RIER, PHIxL, PHIxR, PHIyL, PHIyR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double PI, H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R, detA;
  double L_u, L_p, L_rho;

  double u_t_mat, p_t_mat;

  double shk_spd, shk_u_s, shk_u_L, shk_u_R, zeta = (gamma-1.0)/(gamma+1.0), zts = zeta*zeta;
  double C, rho_x;

  double g_rho, g_u, g_p, f;

  double speed_L, speed_R;


  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);
  cricit= u_R - u_L - 2.0*(c_R+c_L)/(gamma-1.0);

/*
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum(wave_speed, D, U, 0.0, rho_L, rho_R, d_rho_L, d_rho_R, u_L, u_R, d_u_L, d_u_R, v_L, v_R, d_v_L, d_v_R, p_L, p_R, d_p_L, d_p_R, gamma, eps);
    return;
  }
*/
  double g0, g1, g2, g3, g4;


  dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));
//=========acoustic case==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;
    u_star = 0.5*(u_R+u_L);

    if(u_star > lambda)
    {
      U[0] = rho_L;
      U[1] = u_L;
      U[2] = v_L;
      U[3] = p_L;
      c_star = c_L;
      c_square = c_star*c_star;
      if(wave_speed[0] > lambda)
      {
	D[0] = -U[1]*d_rho_L - U[2]*t_rho_L - U[0]*(d_u_L+t_v_L);
	D[1] = -U[1]*  d_u_L - U[2]*  t_u_L - d_p_L/U[0];
	D[2] = -U[1]*  d_v_L - U[2]*  t_v_L - t_p_L/U[0];
	D[3] = -U[1]*  d_p_L - U[2]*  t_p_L - U[0]*(d_u_L+t_v_L)*c_square;
      }
      else
      {
	PHIxL = d_u_L + d_p_L/(U[0]*c_star);
	PHIxR = d_u_R - d_p_R/(U[0]*c_star);
	PHIyL = t_u_L + t_p_L/(U[0]*c_star);
	PHIyR = t_u_R - t_p_R/(U[0]*c_star);
	D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
	D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
	D[2] = -U[1]*d_v_L - U[2]*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
	D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L) + U[2]*(t_p_L - c_square*t_rho_L))/c_square;
      }
      D[0] = D[0] + lambda * d_rho_L;
      D[1] = D[1] + lambda * d_u_L;
      D[2] = D[2] + lambda * d_v_L;
      D[3] = D[3] + lambda * d_p_L;
    }
    else
    {
      U[0] = rho_R;
      U[1] = u_R;
      U[2] = v_R;
      U[3] = p_R;
      c_star = c_R;
      c_square = c_star*c_star;
      if(wave_speed[1] < lambda)
      {
	D[0] = -U[1]*d_rho_R - U[2]*t_rho_R - U[0]*(d_u_R+t_v_R);
	D[1] = -U[1]*  d_u_R - U[2]*  t_u_R - d_p_R/U[0];
	D[2] = -U[1]*  d_v_R - U[2]*  t_v_R - t_p_R/U[0];
	D[3] = -U[1]*  d_p_R - U[2]*  t_p_R - U[0]*(d_u_R+t_v_R)*c_square;
      }
      else
      {
	PHIxL = d_u_L + d_p_L/(U[0]*c_star);
	PHIxR = d_u_R - d_p_R/(U[0]*c_star);
	PHIyL = t_u_L + t_p_L/(U[0]*c_star);
	PHIyR = t_u_R - t_p_R/(U[0]*c_star);
	D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
	D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
	D[2] = -U[1]*d_v_R - U[2]*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
	D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R) + U[2]*(t_p_R - c_square*t_rho_R))/c_square;
      }
      D[0] = D[0] + lambda * d_rho_R;
      D[1] = D[1] + lambda * d_u_R;
      D[2] = D[2] + lambda * d_v_R;
      D[3] = D[3] + lambda * d_p_R;
    }
    return;
  }

//=========non-acoustic case==========
  Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, 500);

  if(CRW[0]){
    rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
      c_star_L =   c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
       speed_L =   u_L - c_L;
  }  else{
    rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
      c_star_L = sqrt(gamma * p_star / rho_star_L);
       speed_L = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
  }
  if(CRW[1]){
    rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
      c_star_R =   c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
       speed_R =   u_R + c_R;
  }  else{
    rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
      c_star_R = sqrt(gamma * p_star / rho_star_R);
       speed_R = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
  }
  wave_speed[0] = speed_L;
  wave_speed[1] = speed_R;

  //------trivial case------
  if(speed_L > lambda) //the direction is on the left side of all the three waves
  {
    U[0] = rho_L;
    U[1] =   u_L;
    U[2] =   v_L;
    U[3] =   p_L;
    c_star = c_L;
    c_square = c_star*c_star;

    D[0] = -U[1]*d_rho_L - U[2]*t_rho_L - U[0]*(d_u_L+t_v_L);
    D[1] = -U[1]*  d_u_L - U[2]*  t_u_L - d_p_L/U[0];
    D[2] = -U[1]*  d_v_L - U[2]*  t_v_L - t_p_L/U[0];
    D[3] = -U[1]*  d_p_L - U[2]*  t_p_L - U[0]*(d_u_L+t_v_L)*c_square;

    D[0] = D[0] + lambda * d_rho_L;
    D[1] = D[1] + lambda * d_u_L;
    D[2] = D[2] + lambda * d_v_L;
    D[3] = D[3] + lambda * d_p_L;
  }
  else if(speed_R < lambda) //the direction is on the right side of all the three waves
  {
    U[0] = rho_R;
    U[1] = u_R;
    U[2] = v_R;
    U[3] = p_R;
    c_star = c_R;
    c_square = c_star*c_star;

    D[0] = -U[1]*d_rho_R - U[2]*t_rho_R - U[0]*(d_u_R+t_v_R);
    D[1] = -U[1]*  d_u_R - U[2]*  t_u_R - d_p_R/U[0];
    D[2] = -U[1]*  d_v_R - U[2]*  t_v_R - t_p_R/U[0];
    D[3] = -U[1]*  d_p_R - U[2]*  t_p_R - U[0]*(d_u_R+t_v_R)*c_square;

    D[0] = D[0] + lambda * d_rho_R;
    D[1] = D[1] + lambda * d_u_R;
    D[2] = D[2] + lambda * d_v_R;
    D[3] = D[3] + lambda * d_p_R;
  }
  else//----non-trivial case----
  {
    if((CRW[0]) && ((u_star-c_star_L) > lambda)) // the direction is in a 1-CRW
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      c_star = U[1] - lambda;
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/c_square;
      U[2] = v_L;

      c_frac = c_star/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -U[1]*d_v_L*U[0]/rho_L;

      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);
      g1 = -0.5*U[2]*(PHIyL + PHIyR);
      g3 = -0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
      g2 = -U[2]*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
      g0 = (g3 + U[2]*(t_p_L - c_square*t_rho_L))/c_square;
    }
    else if((CRW[1]) && ((u_star+c_star_R) < lambda)) // the direction is in a 3-CRW
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      c_star = lambda-U[1];
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/c_square;
      U[2] = v_R;

      c_frac = c_star/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -U[1]*d_v_R*U[0]/rho_R;

      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);
      g1 = -0.5*U[2]*(PHIyL + PHIyR);
      g3 = -0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
      g2 = -U[2]*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
      g0 = (g3 + U[2]*(t_p_R - c_square*t_rho_R))/c_square;
    }
    else//--non-sonic case--
    {
      //determine a_L, b_L and d_L
      if(CRW[0]) //the 1-wave is a CRW
      {
	a_L = 1.0;
        b_L = 1.0 / rho_star_L / c_star_L;
	c_frac = c_star_L/c_L;
	TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
	d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;
	d_L = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
	d_L = d_L * TdS;
	d_L = d_L - c_L*pow(c_frac, 0.5/zeta) * d_Psi;
      }
      else //the 1-wave is a shock
      {
	shk_u_s = -sqrt(0.5*((gamma+1.0)*p_L   +(gamma-1.0)*p_star)/rho_star_L);
	shk_u_L = -sqrt(0.5*((gamma+1.0)*p_star+(gamma-1.0)*p_L   )/rho_L);

	VAR = sqrt((1-zeta)/(rho_L*(p_star+zeta*p_L)));
	H1 =  0.5* VAR * (p_star+(1.0+2.0*zeta)*p_L)/(p_star+zeta*p_L);
	H2 = -0.5*VAR * ((2.0+zeta)*p_star + zeta*p_L)/(p_star+zeta*p_L);
	H3 = -0.5*VAR * (p_star-p_L) / rho_L;

	L_p = -1.0/rho_L - shk_u_L*H2;
	L_u = shk_u_L + rho_L*(c_L*c_L*H2 + H3);
	L_rho = -shk_u_L * H3;

	a_L = 1.0 - rho_star_L* shk_u_s * H1;
	b_L = -shk_u_s/(rho_star_L*c_star_L*c_star_L)+ H1;
	d_L = L_rho*d_rho_L + L_u*d_u_L + L_p*d_p_L;
      }
      //determine a_R, b_R and d_R
      if(CRW[1]) //the 3-wave is a CRW
      {
	a_R = 1.0;
        b_R = -1.0 / rho_star_R / c_star_R;
	c_frac = c_star_R/c_R;
	TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
	d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;
	d_R = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
	d_R = d_R * TdS;
	d_R = d_R + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      }
      else //the 3-wave is a shock
      {
	shk_u_s = sqrt(0.5*((gamma+1.0)*p_R   + (gamma-1.0)*p_star)/rho_star_R);
	shk_u_R = sqrt(0.5*((gamma+1.0)*p_star+ (gamma-1.0)*p_R   )/rho_R);

	VAR  = sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R)));

	H1 = 0.5* VAR * (p_star+(1+2.0*zeta)*p_R)/(p_star+zeta*p_R);
	H2 = -0.5*VAR * ((2.0+zeta)*p_star+zeta*p_R)/(p_star+zeta*p_R);
	H3 = -0.5*(p_star-p_R)* VAR /rho_R;

	L_p = -1.0/rho_R + shk_u_R*H2;
	L_u = shk_u_R - rho_R*(c_R*c_R*H2 + H3);
	L_rho = shk_u_R * H3;

	a_R = 1.0 +rho_star_R* shk_u_s * H1;
	b_R = -(shk_u_s/(rho_star_R*c_star_R*c_star_R) + H1);
	d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R;
      }

      detA = a_L*b_R - b_L*a_R;
      u_t_mat = (b_R*d_L - b_L*d_R)/detA;
      p_t_mat = (a_L*d_R - a_R*d_L)/detA;

      if(u_star < lambda) //the direction is between the contact discontinuety and the 3-wave
      {
	U[0] = rho_star_R;
	U[1] =   u_star;
	U[2] =   v_R;
	U[3] =   p_star;
	c_star = c_star_R;
	c_square = c_star*c_star;

	//already total D!
        D[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_R/c_square;
        D[3] = p_t_mat + rho_star_R*(u_star-lambda) * u_t_mat;

	if(CRW[1]) //the 3-wave is a CRW
	{
	  //already total D!
	  D[0] = rho_star_R*(u_star-lambda)*pow(c_star_R/c_R, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
	  D[0] = (D[0] + D[3]) / c_square;

	  D[2] = -U[1]*d_v_R*U[0]/rho_R;
	  D[2] = D[2] + lambda*d_v_R;
	}
	else //the 3-wave is a shock
	{
	  shk_u_s = sqrt(0.5*((gamma+1.0)*p_R   + (gamma-1.0)*p_star)/rho_star_R);
	  shk_u_R = sqrt(0.5*((gamma+1.0)*p_star+ (gamma-1.0)*p_R   )/rho_R);

	  VAR = p_R + zeta*p_star;
	  H1 = rho_R * p_R    * (1.0 - zts) / (VAR*VAR);
	  H2 = rho_R * p_star * (zts - 1.0) / (VAR*VAR);
	  H3 = (p_star + zeta*p_R)/VAR;

	  L_rho = shk_u_R * H3 * d_rho_R;
	  L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
	  L_p = H2 * shk_u_R * d_p_R;

	  D[0] = ((u_star+shk_u_s)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*shk_u_s*H1*u_t_mat;
	  D[0] = (D[0] - u_star*(L_p+L_rho+L_u)) / shk_u_s;


	  f = shk_u_R*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
	  rho_x = (f + H1*(p_t_mat - rho_star_R*shk_u_s*u_t_mat) - D[0]) / (shk_u_R+u_R);
	  D[0] = D[0] + lambda*rho_x;

	  D[2] = U[1] * shk_u_R * d_v_R / shk_u_s;
	  D[2] = -D[2] + lambda*d_v_R;
	}

	PHIyL = t_u_L + t_p_L/(U[0]*c_star);
	PHIyR = t_u_R - t_p_R/(U[0]*c_star);
        g1 = -0.5*U[2]*(PHIyL + PHIyR);
        g3 = -0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
        g2 = -U[2]*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
        g0 = (g3 + U[2]*(t_p_R - c_square*t_rho_R))/c_square;
      }
      else //the direction is between the 1-wave and the contact discontinuety
      {
	U[0] = rho_star_L;
	U[1] =   u_star;
	U[2] =   v_L;
	U[3] =   p_star;
	c_star = c_star_L;
	c_square = c_star*c_star;

	//already total D!
        D[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_L/c_square;
        D[3] = p_t_mat + rho_star_L*(u_star-lambda) * u_t_mat;
	if(CRW[0]) //the 1-wave is a CRW
	{
	  //already total D!
	  D[0] = rho_star_L*(u_star-lambda)*pow(c_star_L/c_L, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
	  D[0] = (D[0] + D[3]) / c_square;

	  D[2] = -U[1]*d_v_L*U[0]/rho_L;
	  D[2] = D[2] + lambda*d_v_L;
	}
	else //the 1-wave is a shock
	{
	  shk_u_s = -sqrt(0.5*((gamma+1.0)*p_L   +(gamma-1.0)*p_star)/rho_star_L);
	  shk_u_L = -sqrt(0.5*((gamma+1.0)*p_star+(gamma-1.0)*p_L   )/rho_L);

	  VAR = p_L + zeta*p_star;
	  H1 = rho_L * p_L    * (1.0 - zts) / (VAR*VAR);
	  H2 = rho_L * p_star * (zts - 1.0) / (VAR*VAR);
	  H3 = (p_star + zeta*p_L)/VAR;

	  L_rho = shk_u_L * H3 * d_rho_L;
	  L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
	  L_p = H2 * shk_u_L * d_p_L;

	  D[0] = ((u_star+shk_u_s)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*shk_u_s*H1*u_t_mat;
	  D[0] = (D[0] - u_star*(L_p+L_rho+L_u))/ shk_u_s;

	  f = shk_u_L*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
	  rho_x = (f + H1*(p_t_mat - rho_star_L*shk_u_s*u_t_mat) - D[0]) / (shk_u_L+u_L);
	  D[0] = D[0] + lambda*rho_x;

	  D[2] = U[1] * shk_u_L * d_v_L / shk_u_s;
	  D[2] = -D[2] + lambda*d_v_L;
	}

	PHIyL = t_u_L + t_p_L/(U[0]*c_star);
	PHIyR = t_u_R - t_p_R/(U[0]*c_star);
	g1 = -0.5*U[2]*(PHIyL + PHIyR);
        g3 = -0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
        g2 = -U[2]*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
	g0 = (g3 + U[2]*(t_p_L - c_square*t_rho_L))/c_square;
      }
    //--end of non-sonic case--
    }

    D[0] = D[0] + g0;
    D[1] = D[1] + g1;
    D[2] = D[2] + g2;
    D[3] = D[3] + g3;
  //----end of non-trivial case----
  }
}




void vacuum
(double *wave_speed, double *D, double *U, double lambda,
 double rho_L, double rho_R, double d_rho_L, double d_rho_R,
 double   u_L, double   u_R, double   d_u_L, double   d_u_R,
 double   v_L, double   v_R, double   d_v_L, double   d_v_R,
 double   p_L, double   p_R, double   d_p_L, double   d_p_R,
 double gamma, double eps)
{
  double NU = 1.0/(gamma-1.0), zeta = (gamma-1.0)/(gamma+1.0);

  double TdS, d_Psi, d_Phi, HEAD_L, RIE_SL, HEAD_R, RIE_RR;
  double c_L, c_R, c_star_L, c_star_R, C, c_frac;

  wave_speed[0] = u_L - c_L;
  wave_speed[0] = u_R + c_R;

  if(((rho_R < eps) || (p_R < eps)) && ((rho_L < eps) || (p_L < eps)))
  {
    U[0] = 0.0;
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = 0.0;
    D[0] = 0.0;
    D[1] = 0.0;
    D[2] = 0.0;
    D[3] = 0.0;
  }
  else if((rho_R < eps) && (p_R < eps))//right dry,the left wave is a Rarefaction
  {
    c_L = sqrt(gamma*p_L/rho_L);
    HEAD_L = u_L - c_L;//SLOPES OF CHAR.
    RIE_SL = u_L + 2.0*NU*c_L;//Riemann Invariant

    if(lambda < HEAD_L)//x lie to the left of fan-region of
    {
      D[0] = -d_rho_L*u_L - rho_L*d_u_L;
      D[1] = -u_L*d_u_L - d_p_L/rho_L;
      D[2] = -u_L*d_v_L;
      D[3] = -(rho_L*c_L*c_L*d_u_L + u_L*d_p_L);

      D[0] = D[0] + lambda * d_rho_L;
      D[1] = D[1] + lambda * d_u_L;
      D[2] = D[2] + lambda * d_v_L;
      D[3] = D[3] + lambda * d_p_L;

      U[0] = rho_L;
      U[1] =   u_L;
      U[2] =   v_L;
      U[3] =   p_L;
    }
    else if(lambda < RIE_SL)//x lie in the  fan-region: SONIC CASE
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      C = U[1] - lambda;
      U[3] = pow(C/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/C/C;

      c_frac = C/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[3]) /C/C;

      U[2] = v_L;
      D[2] = -U[1]*d_v_L*U[0]/rho_L;
    }	 
    else//x lie at the right of fan-region 
    {
      U[0] = 0.0;
      U[1] = 0.0;
      U[2] = 0.0;
      U[3] = 0.0;
      D[0] = 0.0;
      D[1] = 0.0;
      D[2] = 0.0;
      D[3] = 0.0;
    }
  }
  else if(rho_L < eps || p_L < eps)//left dry,the right wave is a Rarefaction
  {
    c_R = sqrt(gamma*p_R/rho_R);
    HEAD_R = u_R + c_R;
    RIE_RR = u_R - 2.0*NU*c_R;//Riemann Invariant

    if(HEAD_R < lambda)
    {
      D[0] = -d_rho_R*u_R - rho_R*d_u_R;
      D[1] = -u_R*d_u_R - d_p_R/rho_R;
      D[2] = -u_R*d_v_R;
      D[3] = -(rho_R*c_R*c_R*d_u_R + u_R*d_p_R);

      D[0] = D[0] + lambda * d_rho_R;
      D[1] = D[1] + lambda * d_u_R;
      D[2] = D[2] + lambda * d_v_R;
      D[3] = D[3] + lambda * d_p_R;

      U[0] = rho_R;
      U[1] =   u_R;
      U[2] =   v_R;
      U[3] =   p_R;
    }
    else if(RIE_RR < lambda)//sonic case
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      C = lambda-U[1];
      U[3] = pow(C/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/C/C;

      c_frac = C/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[3]) /C/C;

      U[2] = v_R;
      D[2] = -U[1]*d_v_L*U[0]/rho_R;
    }
    else
    {
      U[0] = 0.0;
      U[1] = 0.0;
      U[2] = 0.0;
      U[3] = 0.0;
      D[0] = 0.0;
      D[1] = 0.0;
      D[2] = 0.0;
      D[3] = 0.0;
    }
  }
  else//middle vacuum
  {
    c_R = sqrt(gamma*p_R/rho_R);
    c_L = sqrt(gamma*p_L/rho_L);
    HEAD_L = u_L - c_L;
    RIE_SL = u_L + 2.0*NU*c_L;
    HEAD_R = u_R + c_R;
    RIE_RR = u_R - 2.0*NU*c_R;

    if(lambda < HEAD_L)
    {
      D[0] = -d_rho_L*u_L - rho_L*d_u_L;
      D[1] = -u_L*d_u_L - d_p_L/rho_L;
      D[2] = -u_L*d_v_L;
      D[3] = -(rho_L*c_L*c_L*d_u_L + u_L*d_p_L);

      D[0] = D[0] + lambda * d_rho_L;
      D[1] = D[1] + lambda * d_u_L;
      D[2] = D[2] + lambda * d_v_L;
      D[3] = D[3] + lambda * d_p_L;

      U[0] = rho_L;
      U[1] =   u_L;
      U[2] =   v_L;
      U[3] =   p_L;	
    }
    else if(lambda < RIE_SL)
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      C = U[1] - lambda;
      U[3] = pow(C/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/C/C;

      c_frac = C/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[3]) /C/C;

      U[2] = v_L;
      D[2] = -U[1]*d_v_L*U[0]/rho_L;
    }
    else if(HEAD_R < lambda)
    {
      D[0] = -d_rho_R*u_R - rho_R*d_u_R;
      D[1] = -u_R*d_u_R - d_p_R/rho_R;
      D[2] = -u_R*d_v_R;
      D[3] = -(rho_R*c_R*c_R*d_u_R + u_R*d_p_R);

      D[0] = D[0] + lambda * d_rho_R;
      D[1] = D[1] + lambda * d_u_R;
      D[2] = D[2] + lambda * d_v_R;
      D[3] = D[3] + lambda * d_p_R;

      U[0] = rho_R;
      U[1] =   u_R;
      U[2] =   v_R;
      U[3] =   p_R;
    }
    else if(RIE_RR < lambda)
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      C = lambda-U[1];
      U[3] = pow(C/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/C/C;

      c_frac = C/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[3]) /C/C;

      U[2] = v_R;
      D[2] = -U[1]*d_v_L*U[0]/rho_R;
    }
    else
    {
      U[0] = 0.0;
      U[1] = 0.0;
      U[2] = 0.0;
      U[3] = 0.0;
      D[0] = 0.0;
      D[1] = 0.0;
      D[2] = 0.0;
      D[3] = 0.0;
    }
  }
}



void linear_GRP_solver_noT
(double *wave_speed, double *D, double *U, double lambda, double gamma, double eps,
 double rho_L, double u_L, double v_L, double p_L,
 double rho_R, double u_R, double v_R, double p_R,
 double d_rho_L, double d_u_L, double d_v_L, double d_p_L,
 double d_rho_R, double d_u_R, double d_v_R, double d_p_R)
{
  double dist;
  double c_L, c_R, c_frac, cricit, c_star, c_square;
  int CRW[2];
  double d_Phi, d_Psi, TdS, VAR, RIEL, RIER, PHIxL, PHIxR, PHIyL, PHIyR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double PI, H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R, detA;
  double L_u, L_p, L_rho;

  double u_t_mat, p_t_mat;

  double shk_spd, shk_u_s, shk_u_L, shk_u_R, zeta = (gamma-1.0)/(gamma+1.0), zts = zeta*zeta;
  double C, rho_x;

  double g_rho, g_u, g_p, f;

  double speed_L, speed_R;


  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);
  cricit= u_R - u_L - 2.0*(c_R+c_L)/(gamma-1.0);

/*
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum(wave_speed, D, U, 0.0, rho_L, rho_R, d_rho_L, d_rho_R, u_L, u_R, d_u_L, d_u_R, v_L, v_R, d_v_L, d_v_R, p_L, p_R, d_p_L, d_p_R, gamma, eps);
    return;
  }
*/
  double g0, g1, g2, g3, g4;


  dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));
//=========acoustic case==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;
    u_star = 0.5*(u_R+u_L);

    if(u_star > lambda)
    {
      U[0] = rho_L;
      U[1] = u_L;
      U[2] = v_L;
      U[3] = p_L;
      c_star = c_L;
      c_square = c_star*c_star;
      if(wave_speed[0] > lambda)
      {
	D[0] = -U[1]*d_rho_L - U[0]*d_u_L;
	D[1] = -U[1]*  d_u_L - d_p_L/U[0];
	D[2] = -U[1]*  d_v_L;
	D[3] = -U[1]*  d_p_L - U[0]*d_u_L*c_square;
      }
      else
      {
	PHIxL = d_u_L + d_p_L/(U[0]*c_star);
	PHIxR = d_u_R - d_p_R/(U[0]*c_star);
	D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR));
	D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR))*U[0]*c_star;
	D[2] = -U[1]*d_v_L;
	D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L))/c_square;
      }
      D[0] = D[0] + lambda * d_rho_L;
      D[1] = D[1] + lambda * d_u_L;
      D[2] = D[2] + lambda * d_v_L;
      D[3] = D[3] + lambda * d_p_L;
    }
    else
    {
      U[0] = rho_R;
      U[1] = u_R;
      U[2] = v_R;
      U[3] = p_R;
      c_star = c_R;
      c_square = c_star*c_star;
      if(wave_speed[1] < lambda)
      {
	D[0] = -U[1]*d_rho_R - U[0]*d_u_R;
	D[1] = -U[1]*  d_u_R - d_p_R/U[0];
	D[2] = -U[1]*  d_v_R;
	D[3] = -U[1]*  d_p_R - U[0]*d_u_R*c_square;
      }
      else
      {
	PHIxL = d_u_L + d_p_L/(U[0]*c_star);
	PHIxR = d_u_R - d_p_R/(U[0]*c_star);
	D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR));
	D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR))*U[0]*c_star;
	D[2] = -U[1]*d_v_R;
	D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R))/c_square;
      }
      D[0] = D[0] + lambda * d_rho_R;
      D[1] = D[1] + lambda * d_u_R;
      D[2] = D[2] + lambda * d_v_R;
      D[3] = D[3] + lambda * d_p_R;
    }
    return;
  }

//=========non-acoustic case==========
  Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, 500);

  if(CRW[0]){
    rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
      c_star_L =   c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
       speed_L =   u_L - c_L;
  }  else{
    rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
      c_star_L = sqrt(gamma * p_star / rho_star_L);
       speed_L = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
  }
  if(CRW[1]){
    rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
      c_star_R =   c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
       speed_R =   u_R + c_R;
  }  else{
    rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
      c_star_R = sqrt(gamma * p_star / rho_star_R);
       speed_R = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
  }
  wave_speed[0] = speed_L;
  wave_speed[1] = speed_R;

  //------trivial case------
  if(speed_L > lambda) //the direction is on the left side of all the three waves
  {
    U[0] = rho_L;
    U[1] =   u_L;
    U[2] =   v_L;
    U[3] =   p_L;
    c_star = c_L;
    c_square = c_star*c_star;

    D[0] = -U[1]*d_rho_L - U[0]*d_u_L;
    D[1] = -U[1]*  d_u_L - d_p_L/U[0];
    D[2] = -U[1]*  d_v_L;
    D[3] = -U[1]*  d_p_L - U[0]*d_u_L*c_square;

    D[0] = D[0] + lambda * d_rho_L;
    D[1] = D[1] + lambda * d_u_L;
    D[2] = D[2] + lambda * d_v_L;
    D[3] = D[3] + lambda * d_p_L;
  }
  else if(speed_R < lambda) //the direction is on the right side of all the three waves
  {
    U[0] = rho_R;
    U[1] = u_R;
    U[2] = v_R;
    U[3] = p_R;
    c_star = c_R;
    c_square = c_star*c_star;

    D[0] = -U[1]*d_rho_R - U[0]*d_u_R;
    D[1] = -U[1]*  d_u_R - d_p_R/U[0];
    D[2] = -U[1]*  d_v_R;
    D[3] = -U[1]*  d_p_R - U[0]*d_u_R*c_square;

    D[0] = D[0] + lambda * d_rho_R;
    D[1] = D[1] + lambda * d_u_R;
    D[2] = D[2] + lambda * d_v_R;
    D[3] = D[3] + lambda * d_p_R;
  }
  else//----non-trivial case----
  {
    if((CRW[0]) && ((u_star-c_star_L) > lambda)) // the direction is in a 1-CRW
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      c_star = U[1] - lambda;
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/c_square;
      U[2] = v_L;

      c_frac = c_star/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -U[1]*d_v_L*U[0]/rho_L;
    }
    else if((CRW[1]) && ((u_star+c_star_R) < lambda)) // the direction is in a 3-CRW
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      c_star = lambda-U[1];
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/c_square;
      U[2] = v_R;

      c_frac = c_star/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      D[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      D[1] = D[1] * TdS;
      D[1] = D[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      D[3] = U[0]*(U[1]-lambda)*D[1];

      D[0] = U[0]*(U[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -U[1]*d_v_R*U[0]/rho_R;
    }
    else//--non-sonic case--
    {
      //determine a_L, b_L and d_L
      if(CRW[0]) //the 1-wave is a CRW
      {
	a_L = 1.0;
        b_L = 1.0 / rho_star_L / c_star_L;
	c_frac = c_star_L/c_L;
	TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
	d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;
	d_L = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
	d_L = d_L * TdS;
	d_L = d_L - c_L*pow(c_frac, 0.5/zeta) * d_Psi;
      }
      else //the 1-wave is a shock
      {
	shk_u_s = -sqrt(0.5*((gamma+1.0)*p_L   +(gamma-1.0)*p_star)/rho_star_L);
	shk_u_L = -sqrt(0.5*((gamma+1.0)*p_star+(gamma-1.0)*p_L   )/rho_L);

	VAR = sqrt((1-zeta)/(rho_L*(p_star+zeta*p_L)));
	H1 =  0.5* VAR * (p_star+(1.0+2.0*zeta)*p_L)/(p_star+zeta*p_L);
	H2 = -0.5*VAR * ((2.0+zeta)*p_star + zeta*p_L)/(p_star+zeta*p_L);
	H3 = -0.5*VAR * (p_star-p_L) / rho_L;

	L_p = -1.0/rho_L - shk_u_L*H2;
	L_u = shk_u_L + rho_L*(c_L*c_L*H2 + H3);
	L_rho = -shk_u_L * H3;

	a_L = 1.0 - rho_star_L* shk_u_s * H1;
	b_L = -shk_u_s/(rho_star_L*c_star_L*c_star_L)+ H1;
	d_L = L_rho*d_rho_L + L_u*d_u_L + L_p*d_p_L;
      }
      //determine a_R, b_R and d_R
      if(CRW[1]) //the 3-wave is a CRW
      {
	a_R = 1.0;
        b_R = -1.0 / rho_star_R / c_star_R;
	c_frac = c_star_R/c_R;
	TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
	d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;
	d_R = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
	d_R = d_R * TdS;
	d_R = d_R + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      }
      else //the 3-wave is a shock
      {
	shk_u_s = sqrt(0.5*((gamma+1.0)*p_R   + (gamma-1.0)*p_star)/rho_star_R);
	shk_u_R = sqrt(0.5*((gamma+1.0)*p_star+ (gamma-1.0)*p_R   )/rho_R);

	VAR  = sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R)));

	H1 = 0.5* VAR * (p_star+(1+2.0*zeta)*p_R)/(p_star+zeta*p_R);
	H2 = -0.5*VAR * ((2.0+zeta)*p_star+zeta*p_R)/(p_star+zeta*p_R);
	H3 = -0.5*(p_star-p_R)* VAR /rho_R;

	L_p = -1.0/rho_R + shk_u_R*H2;
	L_u = shk_u_R - rho_R*(c_R*c_R*H2 + H3);
	L_rho = shk_u_R * H3;

	a_R = 1.0 +rho_star_R* shk_u_s * H1;
	b_R = -(shk_u_s/(rho_star_R*c_star_R*c_star_R) + H1);
	d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R;
      }

      detA = a_L*b_R - b_L*a_R;
      u_t_mat = (b_R*d_L - b_L*d_R)/detA;
      p_t_mat = (a_L*d_R - a_R*d_L)/detA;

      if(u_star < lambda) //the direction is between the contact discontinuety and the 3-wave
      {
	U[0] = rho_star_R;
	U[1] =   u_star;
	U[2] =   v_R;
	U[3] =   p_star;
	c_star = c_star_R;
	c_square = c_star*c_star;

	//already total D!
        D[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_R/c_square;
        D[3] = p_t_mat + rho_star_R*(u_star-lambda) * u_t_mat;

	if(CRW[1]) //the 3-wave is a CRW
	{
	  //already total D!
	  D[0] = rho_star_R*(u_star-lambda)*pow(c_star_R/c_R, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
	  D[0] = (D[0] + D[3]) / c_square;

	  D[2] = -U[1]*d_v_R*U[0]/rho_R;
	  D[2] = D[2] + lambda*d_v_R;
	}
	else //the 3-wave is a shock
	{
	  shk_u_s = sqrt(0.5*((gamma+1.0)*p_R   + (gamma-1.0)*p_star)/rho_star_R);
	  shk_u_R = sqrt(0.5*((gamma+1.0)*p_star+ (gamma-1.0)*p_R   )/rho_R);

	  VAR = p_R + zeta*p_star;
	  H1 = rho_R * p_R    * (1.0 - zts) / (VAR*VAR);
	  H2 = rho_R * p_star * (zts - 1.0) / (VAR*VAR);
	  H3 = (p_star + zeta*p_R)/VAR;

	  L_rho = shk_u_R * H3 * d_rho_R;
	  L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
	  L_p = H2 * shk_u_R * d_p_R;

	  D[0] = ((u_star+shk_u_s)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*shk_u_s*H1*u_t_mat;
	  D[0] = (D[0] - u_star*(L_p+L_rho+L_u)) / shk_u_s;


	  f = shk_u_R*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
	  rho_x = (f + H1*(p_t_mat - rho_star_R*shk_u_s*u_t_mat) - D[0]) / (shk_u_R+u_R);
	  D[0] = D[0] + lambda*rho_x;

	  D[2] = U[1] * shk_u_R * d_v_R / shk_u_s;
	  D[2] = -D[2] + lambda*d_v_R;
	}
      }
      else //the direction is between the 1-wave and the contact discontinuety
      {
	U[0] = rho_star_L;
	U[1] =   u_star;
	U[2] =   v_L;
	U[3] =   p_star;
	c_star = c_star_L;
	c_square = c_star*c_star;

	//already total D!
        D[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_L/c_square;
        D[3] = p_t_mat + rho_star_L*(u_star-lambda) * u_t_mat;
	if(CRW[0]) //the 1-wave is a CRW
	{
	  //already total D!
	  D[0] = rho_star_L*(u_star-lambda)*pow(c_star_L/c_L, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
	  D[0] = (D[0] + D[3]) / c_square;

	  D[2] = -U[1]*d_v_L*U[0]/rho_L;
	  D[2] = D[2] + lambda*d_v_L;
	}
	else //the 1-wave is a shock
	{
	  shk_u_s = -sqrt(0.5*((gamma+1.0)*p_L   +(gamma-1.0)*p_star)/rho_star_L);
	  shk_u_L = -sqrt(0.5*((gamma+1.0)*p_star+(gamma-1.0)*p_L   )/rho_L);

	  VAR = p_L + zeta*p_star;
	  H1 = rho_L * p_L    * (1.0 - zts) / (VAR*VAR);
	  H2 = rho_L * p_star * (zts - 1.0) / (VAR*VAR);
	  H3 = (p_star + zeta*p_L)/VAR;

	  L_rho = shk_u_L * H3 * d_rho_L;
	  L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
	  L_p = H2 * shk_u_L * d_p_L;

	  D[0] = ((u_star+shk_u_s)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*shk_u_s*H1*u_t_mat;
	  D[0] = (D[0] - u_star*(L_p+L_rho+L_u))/ shk_u_s;

	  f = shk_u_L*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
	  rho_x = (f + H1*(p_t_mat - rho_star_L*shk_u_s*u_t_mat) - D[0]) / (shk_u_L+u_L);
	  D[0] = D[0] + lambda*rho_x;

	  D[2] = U[1] * shk_u_L * d_v_L / shk_u_s;
	  D[2] = -D[2] + lambda*d_v_L;
	}
      }
    //--end of non-sonic case--
    }
  //----end of non-trivial case----
  }
}
