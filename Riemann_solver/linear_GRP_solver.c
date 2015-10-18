#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"



void linear_GRP_solver
(double *wave_speed, double *direvative, double *source, double lambda,
 double rho_L, double rho_R, double d_rho_L, double d_rho_R,
 double   u_L, double   u_R, double   d_u_L, double   d_u_R,
 double   v_L, double   v_R, double   d_v_L, double   d_v_R,
 double   p_L, double   p_R, double   d_p_L, double   d_p_R,
 double gamma, double eps)
{
  double dist;
  double c_L, c_R, c_frac, cricit;
  int CRW[2];
  double d_Phi, d_Psi, TdS, VAR, RIEL, RIER;
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


  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum(wave_speed, direvative, source, 0.0, rho_L, rho_R, d_rho_L, d_rho_R, u_L, u_R, d_u_L, d_u_R, v_L, v_R, d_v_L, d_v_R, p_L, p_R, d_p_L, d_p_R, gamma, eps);
    return;
  }


  dist = (rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  //dist = sqrt((u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R));
//=========acoustic case==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;

    rho_star_L = rho_L;
    rho_star_R = rho_R;
      c_star_L =   c_L;
      c_star_R =   c_R;
      u_star = 0.5*(u_R+u_L);
      p_star = 0.5*(p_R+p_L);

    RIEL = d_u_L + d_p_L/(rho_L*c_L);
    RIER = d_u_R - d_p_R/(rho_R*c_R);
    direvative[1] = -0.5 * ((u_L+c_L)*RIEL + (u_R-c_R)*RIER);
    direvative[3] = -0.5 * rho_L*c_L*((u_L+c_L)*RIEL - (u_R-c_R)*RIER);

    if(u_star > lambda)
    {
      source[0] = rho_star_L;
      source[1] =   u_star;
      source[2] =   v_L;
      source[3] =   p_star;
      direvative[0] = (direvative[1] + u_star*(d_p_L - c_L*c_L*d_rho_L)) /c_L/c_L;
      direvative[2] = -u_L*d_v_L;

      direvative[0] = direvative[0] + lambda * d_rho_L;
      direvative[1] = direvative[1] + lambda * d_u_L;
      direvative[2] = direvative[2] + lambda * d_v_L;
      direvative[3] = direvative[3] + lambda * d_p_L;
    }
    else
    {
      source[0] = rho_star_R;
      source[1] =   u_star;
      source[2] =   v_R;
      source[3] =   p_star;
      direvative[0] = (direvative[1] + u_star*(d_p_R - c_R*c_R*d_rho_R)) /c_R/c_R;
      direvative[2] = -u_R*d_v_R;

      direvative[0] = direvative[0] + lambda * d_rho_R;
      direvative[1] = direvative[1] + lambda * d_u_R;
      direvative[2] = direvative[2] + lambda * d_v_R;
      direvative[3] = direvative[3] + lambda * d_p_R;
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
    direvative[0] = -d_rho_L*u_L - rho_L*d_u_L;
    direvative[1] = -u_L*d_u_L - d_p_L/rho_L;
    direvative[2] = -u_L*d_v_L;
    direvative[3] = -(rho_L*c_L*c_L*d_u_L + u_L*d_p_L);

    direvative[0] = direvative[0] + lambda * d_rho_L;
    direvative[1] = direvative[1] + lambda * d_u_L;
    direvative[2] = direvative[2] + lambda * d_v_L;
    direvative[3] = direvative[3] + lambda * d_p_L;

    source[0] = rho_L;
    source[1] =   u_L;
    source[2] =   v_L;
    source[3] =   p_L;
  }
  else if(speed_R < lambda) //the direction is on the right side of all the three waves
  {
    direvative[0] = -d_rho_R*u_R - rho_R*d_u_R;
    direvative[1] = -u_R*d_u_R - d_p_R/rho_R;
    direvative[2] = -u_R*d_v_R;
    direvative[3] = -(rho_R*c_R*c_R*d_u_R + u_R*d_p_R);

    direvative[0] = direvative[0] + lambda * d_rho_R;
    direvative[1] = direvative[1] + lambda * d_u_R;
    direvative[2] = direvative[2] + lambda * d_v_R;
    direvative[3] = direvative[3] + lambda * d_p_R;

    source[0] = rho_R;
    source[1] =   u_R;
    source[2] =   v_R;
    source[3] =   p_R;
  }
  else//----non-trivial case----
  {
    if((CRW[0]) && ((u_star-c_star_L) > lambda)) // the direction is in a 1-CRW
    {
      source[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      C = source[1] - lambda;
      source[3] = pow(C/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_L;
      direvative[2] = -source[1]*d_v_L*source[0]/rho_L;
    }
    else if((CRW[1]) && ((u_star+c_star_R) < lambda)) // the direction is in a 3-CRW
    {
      source[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      C = lambda-source[1];
      source[3] = pow(C/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_R;
      direvative[2] = -source[1]*d_v_R*source[0]/rho_R;
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
	source[0] = rho_star_R;
	source[1] =   u_star;
	source[2] =   v_R;
	source[3] =   p_star;

	//already total direvative!
        direvative[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_R/c_star_R/c_star_R;
        direvative[3] = p_t_mat + rho_star_R*(u_star-lambda) * u_t_mat;

	if(CRW[1]) //the 3-wave is a CRW
	{
	  //already total direvative!
	  direvative[0] = rho_star_R*(u_star-lambda)*pow(c_star_R/c_R, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
	  direvative[0] = (direvative[0] + direvative[3]) / c_star_R/c_star_R;

	  direvative[2] = -source[1]*d_v_R*source[0]/rho_R;
	  direvative[2] = direvative[2] + lambda*d_v_R;
	}
	else //the 3-wave is a shock
	{
	  shk_u_s = sqrt(0.5*((gamma+1.0)*p_R   + (gamma-1.0)*p_star)/rho_star_R);
	  shk_u_R = sqrt(0.5*((gamma+1.0)*p_star+ (gamma-1.0)*p_R   )/rho_R);

	  VAR = p_R + zeta*p_star;
	  H1 = rho_R * p_R    * (1.0 - zts) / VAR/VAR;
	  H2 = rho_R * p_star * (zts - 1.0) / VAR/VAR;
	  H3 = (p_star + zeta*p_R)/VAR;

	  L_rho = shk_u_R * H3 * d_rho_R;
	  L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
	  L_p = H2 * shk_u_R * d_p_R;

	  direvative[0] = ((u_star+shk_u_s)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*shk_u_s*H1*u_t_mat;
	  direvative[0] = (direvative[0] - u_star*(L_p+L_rho+L_u)) / shk_u_s;


	  f = shk_u_R*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
	  rho_x = (f + H1*(p_t_mat - rho_star_R*shk_u_s*u_t_mat) - direvative[0]) / (shk_u_R+u_R);
	  direvative[0] = direvative[0] + lambda*rho_x;

	  direvative[2] = source[1] * shk_u_R * d_v_R / shk_u_s;
	  direvative[2] = -direvative[2] + lambda*d_v_R;
	}
      }
      else //the direction is between the 1-wave and the contact discontinuety
      {
	source[0] = rho_star_L;
	source[1] =   u_star;
	source[2] =   v_L;
	source[3] =   p_star;
	//already total direvative!
        direvative[1] = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_L/c_star_L/c_star_L;
        direvative[3] = p_t_mat + rho_star_L*(u_star-lambda) * u_t_mat;
	if(CRW[0]) //the 1-wave is a CRW
	{
	  //already total direvative!
	  direvative[0] = rho_star_L*(u_star-lambda)*pow(c_star_L/c_L, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
	  direvative[0] = (direvative[0] + direvative[3]) / c_star_L/c_star_L;

	  direvative[2] = -source[1]*d_v_L*source[0]/rho_L;
	  direvative[2] = direvative[2] + lambda*d_v_L;
	}
	else //the 1-wave is a shock
	{
	  shk_u_s = -sqrt(0.5*((gamma+1.0)*p_L   +(gamma-1.0)*p_star)/rho_star_L);
	  shk_u_L = -sqrt(0.5*((gamma+1.0)*p_star+(gamma-1.0)*p_L   )/rho_L);

	  VAR = p_L + zeta*p_star;
	  H1 = rho_L * p_L    * (1.0 - zts) / VAR/VAR;
	  H2 = rho_L * p_star * (zts - 1.0) / VAR/VAR;
	  H3 = (p_star + zeta*p_L)/VAR;

	  L_rho = shk_u_L * H3 * d_rho_L;
	  L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
	  L_p = H2 * shk_u_L * d_p_L;

	  direvative[0] = ((u_star+shk_u_s)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*shk_u_s*H1*u_t_mat;
	  direvative[0] = (direvative[0] - u_star*(L_p+L_rho+L_u))/ shk_u_s;

	  f = shk_u_L*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
	  rho_x = (f + H1*(p_t_mat - rho_star_L*shk_u_s*u_t_mat) - direvative[0]) / (shk_u_L+u_L);
	  direvative[0] = direvative[0] + lambda*rho_x;

	  direvative[2] = source[1] * shk_u_L * d_v_L / shk_u_s;
	  direvative[2] = -direvative[2] + lambda*d_v_L;
	}
      }
    //--end of non-sonic case--
    }
  //----end of non-trivial case----
  }
}


void vacuum
(double *wave_speed, double *direvative, double *source, double lambda,
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
    source[0] = 0.0;
    source[1] = 0.0;
    source[2] = 0.0;
    source[3] = 0.0;
    direvative[0] = 0.0;
    direvative[1] = 0.0;
    direvative[2] = 0.0;
    direvative[3] = 0.0;
  }
  else if((rho_R < eps) && (p_R < eps))//right dry,the left wave is a Rarefaction
  {
    c_L = sqrt(gamma*p_L/rho_L);
    HEAD_L = u_L - c_L;//SLOPES OF CHAR.
    RIE_SL = u_L + 2.0*NU*c_L;//Riemann Invariant

    if(lambda < HEAD_L)//x lie to the left of fan-region of
    {
      direvative[0] = -d_rho_L*u_L - rho_L*d_u_L;
      direvative[1] = -u_L*d_u_L - d_p_L/rho_L;
      direvative[2] = -u_L*d_v_L;
      direvative[3] = -(rho_L*c_L*c_L*d_u_L + u_L*d_p_L);

      direvative[0] = direvative[0] + lambda * d_rho_L;
      direvative[1] = direvative[1] + lambda * d_u_L;
      direvative[2] = direvative[2] + lambda * d_v_L;
      direvative[3] = direvative[3] + lambda * d_p_L;

      source[0] = rho_L;
      source[1] =   u_L;
      source[2] =   v_L;
      source[3] =   p_L;
    }
    else if(lambda < RIE_SL)//x lie in the  fan-region: SONIC CASE
    {
      source[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      C = source[1] - lambda;
      source[3] = pow(C/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_L;
      direvative[2] = -source[1]*d_v_L*source[0]/rho_L;
    }	 
    else//x lie at the right of fan-region 
    {
      source[0] = 0.0;
      source[1] = 0.0;
      source[2] = 0.0;
      source[3] = 0.0;
      direvative[0] = 0.0;
      direvative[1] = 0.0;
      direvative[2] = 0.0;
      direvative[3] = 0.0;
    }
  }
  else if(rho_L < eps || p_L < eps)//left dry,the right wave is a Rarefaction
  {
    c_R = sqrt(gamma*p_R/rho_R);
    HEAD_R = u_R + c_R;
    RIE_RR = u_R - 2.0*NU*c_R;//Riemann Invariant

    if(HEAD_R < lambda)
    {
      direvative[0] = -d_rho_R*u_R - rho_R*d_u_R;
      direvative[1] = -u_R*d_u_R - d_p_R/rho_R;
      direvative[2] = -u_R*d_v_R;
      direvative[3] = -(rho_R*c_R*c_R*d_u_R + u_R*d_p_R);

      direvative[0] = direvative[0] + lambda * d_rho_R;
      direvative[1] = direvative[1] + lambda * d_u_R;
      direvative[2] = direvative[2] + lambda * d_v_R;
      direvative[3] = direvative[3] + lambda * d_p_R;

      source[0] = rho_R;
      source[1] =   u_R;
      source[2] =   v_R;
      source[3] =   p_R;
    }
    else if(RIE_RR < lambda)//sonic case
    {
      source[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      C = lambda-source[1];
      source[3] = pow(C/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_R;
      direvative[2] = -source[1]*d_v_L*source[0]/rho_R;
    }
    else
    {
      source[0] = 0.0;
      source[1] = 0.0;
      source[2] = 0.0;
      source[3] = 0.0;
      direvative[0] = 0.0;
      direvative[1] = 0.0;
      direvative[2] = 0.0;
      direvative[3] = 0.0;
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
      direvative[0] = -d_rho_L*u_L - rho_L*d_u_L;
      direvative[1] = -u_L*d_u_L - d_p_L/rho_L;
      direvative[2] = -u_L*d_v_L;
      direvative[3] = -(rho_L*c_L*c_L*d_u_L + u_L*d_p_L);

      direvative[0] = direvative[0] + lambda * d_rho_L;
      direvative[1] = direvative[1] + lambda * d_u_L;
      direvative[2] = direvative[2] + lambda * d_v_L;
      direvative[3] = direvative[3] + lambda * d_p_L;

      source[0] = rho_L;
      source[1] =   u_L;
      source[2] =   v_L;
      source[3] =   p_L;	
    }
    else if(lambda < RIE_SL)
    {
      source[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      C = source[1] - lambda;
      source[3] = pow(C/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_L;
      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      d_Psi = d_u_L + (gamma*d_p_L/c_L - c_L*d_rho_L)/(gamma-1.0)/rho_L;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta) )/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] - c_L*pow(c_frac, 0.5/zeta)*d_Psi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_L;
      direvative[2] = -source[1]*d_v_L*source[0]/rho_L;
    }
    else if(HEAD_R < lambda)
    {
      direvative[0] = -d_rho_R*u_R - rho_R*d_u_R;
      direvative[1] = -u_R*d_u_R - d_p_R/rho_R;
      direvative[2] = -u_R*d_v_R;
      direvative[3] = -(rho_R*c_R*c_R*d_u_R + u_R*d_p_R);

      direvative[0] = direvative[0] + lambda * d_rho_R;
      direvative[1] = direvative[1] + lambda * d_u_R;
      direvative[2] = direvative[2] + lambda * d_v_R;
      direvative[3] = direvative[3] + lambda * d_p_R;

      source[0] = rho_R;
      source[1] =   u_R;
      source[2] =   v_R;
      source[3] =   p_R;
    }
    else if(RIE_RR < lambda)
    {
      source[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      C = lambda-source[1];
      source[3] = pow(C/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      source[0] = gamma*source[3]/C/C;

      c_frac = C/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      d_Phi = d_u_R - (gamma*d_p_R/c_R - c_R*d_rho_R)/(gamma-1.0)/rho_R;

      direvative[1] = ( (1.0+zeta)*pow(c_frac, 0.5/zeta) + zeta*pow(c_frac, (1.0+zeta)/zeta))/(1.0+2.0*zeta);
      direvative[1] = direvative[1] * TdS;
      direvative[1] = direvative[1] + c_R*pow(c_frac, 0.5/zeta)*d_Phi;
      direvative[3] = source[0]*(source[1]-lambda)*direvative[1];

      direvative[0] = source[0]*(source[1]-lambda)*pow(c_frac, (1.0+zeta)/zeta)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      direvative[0] = (direvative[0] + direvative[3]) /C/C;

      source[2] = v_R;
      direvative[2] = -source[1]*d_v_L*source[0]/rho_R;
    }
    else
    {
      source[0] = 0.0;
      source[1] = 0.0;
      source[2] = 0.0;
      source[3] = 0.0;
      direvative[0] = 0.0;
      direvative[1] = 0.0;
      direvative[2] = 0.0;
      direvative[3] = 0.0;
    }
  }
}
