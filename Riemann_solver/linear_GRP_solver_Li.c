#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"



void linear_GRP_solver_Li(double *UL, double *UR, double *DL, double *DR, 
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   v_L, double   v_R, double   s_v_L, double   s_v_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps)
{
  double NU = 1.0/(gamma-1.0), MU2 = (gamma-1.0)/(gamma+1.0);
  double CL, CR, CRICIT;
  double WL[3], WR[3], SLOPEL[3], SLOPER[3];

  WL[0] = rho_L;
  WL[1] =   u_L;
  WL[2] =   p_L;
  WR[0] = rho_R;
  WR[1] =   u_R;
  WR[2] =   p_R;
  SLOPEL[0] = s_rho_L;
  SLOPEL[1] =   s_u_L;
  SLOPEL[2] =   s_p_L;
  SLOPER[0] = s_rho_R;
  SLOPER[1] =   s_u_R;
  SLOPER[2] =   s_p_R;

  if(WL[0]<eps || WL[2]<eps || WR[0]<eps || WR[2]<eps)
    VACUUM(UL, UR,  DL, DR, gamma, eps, WL, WR, SLOPEL, SLOPER);
  else
  {
    CL=sqrt(gamma*WL[2]/WL[0]);
    CR=sqrt(gamma*WR[2]/WR[0]);
    CRICIT= WR[1]-WL[1]-2.0*NU*(CR+CL);

    if(CRICIT>-eps)
      VACUUM(UL, UR,  DL, DR, gamma, eps, WL, WR, SLOPEL, SLOPER);
    else
      NONVACUUM(UL, UR, DL, DR, gamma, eps, WL, WR, SLOPEL, SLOPER);
  }
}

void VACUUM(double *UL, double *UR, double *DL, double *DR, double GAMMA, double eps, double *WL, double *WR, double *SLOPEL, double *SLOPER)
{
  double NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);

  double A[3][3], K[3], MID[3], PT[3];
  double CL, CR, CSTAR, HEAD_L, RIE_SL, HEAD_R, RIE_RR;

  if((WR[0]<eps || WR[2]<eps) && (WL[0]<eps || WL[2]<eps))
  {
    MID[0] = 0.0;
    MID[1] = 0.0;
    MID[2] = 0.0;
    PT[0] = 0.0;
    PT[1] = 0.0;
    PT[2] = 0.0;
  }
  else if(WR[0]<eps && WR[2]<eps)//right dry,the left wave is a Rarefaction
  {
    CL = sqrt(GAMMA*WL[2]/WL[0]);
    HEAD_L = WL[1] - CL;//SLOPES OF CHAR.
    RIE_SL = WL[1] + 2.0*NU*CL;//Riemann Invariant

    if(0.0<HEAD_L)//x lie to the left of fan-region of
    {
      MID[0] = WL[0];
      MID[1] = WL[1];
      MID[2] = WL[2];
      COEFFICIENT_MATRIX(A,GAMMA,MID);
      MATMUL(PT, A, SLOPEL);//left Rarefaction
      PT[0] = -PT[0];
      PT[1] = -PT[1];
      PT[2] = -PT[2];
    }
    else if(0.0<RIE_SL)//x lie in the  fan-region: SONIC CASE
    {
      CSTAR = MU2* RIE_SL;
      MID[0] = WL[0]*pow((CSTAR/CL),2.0*NU);
      MID[1] = CSTAR;
      MID[2] = WL[2]*pow(MID[0]/WL[0],GAMMA);

      //return d_L which is the time direvative of the velocity
      K[2] = RAREFACTION_LEFT(GAMMA, CSTAR, WL, CL, SLOPEL);
      PT[1] = K[2];
      PT[2] = MID[0] * MID[1] * K[2];
      //return the time direvative of the density
      PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WL, CL, SLOPEL);
    }	 
    else//x lie at the right of fan-region 
    {
      MID[0] = 0.0;
      MID[1] = 0.0;
      MID[2] = 0.0;
      PT[0] = 0.0;
      PT[1] = 0.0;
      PT[2] = 0.0;
    }
  }
  else if(WL[0]<eps || WL[2]<eps)//left dry,the right wave is a Rarefaction
  {
    CR = sqrt(GAMMA*WR[2]/WR[0]);
    HEAD_R = WR[1] + CR;
    RIE_RR = WR[1] - 2.0*NU*CR;//Riemann Invariant

    if(HEAD_R<0.0)
    {
      MID[0] = WR[0];
      MID[1] = WR[1];
      MID[2] = WR[2];
      COEFFICIENT_MATRIX(A, GAMMA, MID);
      MATMUL(PT, A, SLOPER);
      PT[0] = -PT[0];
      PT[1] = -PT[1];
      PT[2] = -PT[2];
    }
    else if(RIE_RR<0.0)//sonic case
    {
      CSTAR = -MU2*RIE_RR;
      MID[1] = -CSTAR;
      MID[0] = WR[0]*pow(CSTAR/CR,2.0*NU);
      MID[2] = WR[2]*pow(MID[0]/WR[0],GAMMA);

      K[2] = RAREFACTION_RIGHT(K[2], CSTAR, WR, CR, SLOPER);
      PT[1] = K[2];
      PT[2] = MID[0]* MID[1] * K[2];
      PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WR, CR, SLOPER);
    }
    else
    {
      MID[0] = 0.0;
      MID[1] = 0.0;
      MID[2] = 0.0;
      PT[0] = 0.0;
      PT[1] = 0.0;
      PT[2] = 0.0;
    }
  }
  else//middle vacuum
  {
    CR = sqrt(GAMMA*WR[2]/WR[0]);
    CL = sqrt(GAMMA*WL[2]/WL[0]);
    HEAD_L = WL[1] - CL;
    RIE_SL = WL[1] + 2.0*NU*CL;
    HEAD_R = WR[1] + CR;
    RIE_RR = WR[1] - 2.0*NU*CR;

    if(0.0<HEAD_L)
    {
      MID[0] = WL[0];
      MID[1] = WL[1];
      MID[2] = WL[2];
      COEFFICIENT_MATRIX(A, GAMMA, MID);
      MATMUL(PT, A, SLOPEL);
      PT[0] = -PT[0];
      PT[1] = -PT[1];
      PT[2] = -PT[2];	
    }
    else if(0.0 < RIE_SL)
    {
      CSTAR = MU2* RIE_SL;
      MID[1] = CSTAR;
      MID[0] = WL[0]*pow(CSTAR/CL,2.0*NU);
      MID[2] = WL[2]*pow(MID[0]/WL[0],GAMMA);

      //return d_L which is the time direvative of the velocity
      K[2] = RAREFACTION_LEFT(GAMMA, CSTAR, WL, CL, SLOPEL);
      PT[1] = K[2];
      PT[2] = MID[0] * MID[1] * K[2];
      //return the time direvative of the density
      PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WL, CL, SLOPEL);
    }
    else if(HEAD_R<0.0)
    {
      MID[0] = WR[0];
      MID[1] = WR[1];
      MID[2] = WR[2];
      COEFFICIENT_MATRIX(A, GAMMA, MID);
      MATMUL(PT, A, SLOPER);
      PT[0] = -PT[0];
      PT[1] = -PT[1];
      PT[2] = -PT[2];
    }
    else if(RIE_RR<0.0)
    {
      CSTAR = -MU2*RIE_RR;
      MID[1] = -CSTAR;
      MID[0] = WR[0]*pow(CSTAR/CR,2.0*NU);
      MID[2] = WR[2]*pow(MID[0]/WR[0],GAMMA);

      K[2] = RAREFACTION_RIGHT(K[2], CSTAR, WR, CR, SLOPER);
      PT[1] = K[2];
      PT[2] = MID[0]* MID[1] * K[2];
      PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WR, CR, SLOPER);
    }
    else
    {
      MID[0] = 0.0;
      MID[1] = 0.0;
      MID[2] = 0.0;
      PT[0] = 0.0;
      PT[1] = 0.0;
      PT[2] = 0.0;
    }
  }
  UL[0] = MID[0];
  UL[1] = MID[1];
  UL[2] = MID[2];
  UR[0] = MID[0];
  UR[1] = MID[1];
  UR[2] = MID[2];
  DL[0] = PT[0];
  DL[1] = PT[1];
  DL[2] = PT[2];
  DR[0] = PT[0];
  DR[1] = PT[1];
  DR[2] = PT[2];
}

void NONVACUUM(double *UL, double *UR, double *DL, double *DR, double GAMMA, double eps, double *WL, double *WR, double *SLOPEL, double *SLOPER)
{
  double NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);

  double A[3][3], B[2], K[3], MID[3], PT[3], WW[3], detA;
  double U_STAR, P_STAR;
  int CRW[2];
  double CL, CR, CSTAR, HEAD_L, RIE_SL, HEAD_R, RIE_RR, TAIL_L, TAIL_R;
  double CSTAR1, CSTAR2, RHO1, RHO2, QL, QR, SL, SR;
  double VAR1, VAR2;
  double G1=(GAMMA-1.0)/(2.0*GAMMA), G2=(GAMMA+1.0)/(2.0*GAMMA);
  double G3=2.0*GAMMA/(GAMMA-1.0), G4=2.0/(GAMMA-1.0);
  double G5=2.0/(GAMMA+1.0), G6=(GAMMA-1.0)/(GAMMA+1.0);
  double G7 = (GAMMA-1.0)/2.0, G8 = GAMMA - 1.0;

  CL = sqrt(GAMMA * WL[2] / WL[0]);
  CR = sqrt(GAMMA * WR[2] / WR[0]);
 	
  RIE_SL = WL[1] + 2.0 * NU * CL;
  RIE_RR = WR[1] - 2.0 * NU * CR;

  VAR1 = (WL[0]-WR[0])*(WL[0]-WR[0]) +  (WL[1]-WR[1])*(WL[1]-WR[1]) + (WL[2]-WR[2])*(WL[2]-WR[2]);

  if(VAR1 < eps)
  {
    ACOUSTIC(PT, WL, CL, WR, CR, SLOPEL, SLOPER);
    UL[0] = WL[0];
    UL[1] = WL[1];
    UL[2] = WL[2];
    UR[0] = WR[0];
    UR[1] = WR[1];
    UR[2] = WR[2];
    DL[0] = PT[0];
    DL[1] = PT[1];
    DL[2] = PT[2];
    DR[0] = PT[0];
    DR[1] = PT[1];
    DR[2] = PT[2];

    return;
  }

  Riemann_solver_exact(&U_STAR, &P_STAR, GAMMA, WL[1], WR[1], WL[2], WR[2], CL, CR, CRW, eps, 50);
  /*
  if(P_STAR > WL[2])
    rho_star_L = rho_L*(P_STAR+zeta*WL[2])/(WL[2]+zeta*P_STAR);
  else
    rho_star_L = rho_L*pow(P_STAR/WL[2],1.0/GAMMA);
  if(P_STAR > WR[2])
    rho_star_R = rho_R*(P_STAR+zeta*WR[2])/(WR[2]+zeta*P_STAR);
  else
    rho_star_R = rho_R*pow(P_STAR/WR[2],1.0/GAMMA);
  c_star_L = sqrt(GAMMA * P_STAR / rho_star_L);
  c_star_R = sqrt(GAMMA * P_STAR / rho_star_R); 
  */

  if(eps < U_STAR)//the t-axe lies to the left of the contact
  {
    if(CRW[0])//left rarefaction
    {
      CSTAR = CL * pow(P_STAR/WL[2] , G1);//G1 = (GAMMA - 1.0)/(2.0*GAMMA)
      HEAD_L = WL[1] - CL;
      TAIL_L = U_STAR - CSTAR;

      if(0.0 < HEAD_L)//left trivial case
      {
	MID[0] = WL[0];
	MID[1] = WL[1];
	MID[2] = WL[2];
	COEFFICIENT_MATRIX(A, GAMMA, MID);
	MATMUL(PT, A, SLOPEL);
	PT[0] = -PT[0];
	PT[1] = -PT[1];
	PT[2] = -PT[2];
      }
      else if(0.0 < TAIL_L)//left sonic case
      {
	CSTAR = MU2*RIE_SL;
	MID[1] = CSTAR;
	MID[0] = WL[0]*pow(CSTAR/CL , 2.0*NU);
	MID[2] = WL[2]*pow(MID[0]/WL[0],GAMMA);
	      
	//return d_L which is the time direvative of the velocity
	K[2] = RAREFACTION_LEFT(GAMMA, CSTAR, WL, CL, SLOPEL);
	PT[1] = K[2];
	PT[2] = MID[0] * MID[1] * K[2];
	//return the time direvative of the density
	PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WL, CL, SLOPEL);
      }
      else// non-sonic case
      {
	MID[2] = P_STAR;
	MID[1] = U_STAR;
	MID[0] = WL[0] * pow(MID[2]/WL[2] , 1.0/GAMMA);
	CSTAR  = sqrt(GAMMA * MID[2]/ MID[0]);
	  
	K[2] = RAREFACTION_LEFT(GAMMA, CSTAR, WL, CL, SLOPEL);
	A[0][0] = 1.0;
	A[0][1] = 1.0/(MID[0]*CSTAR);
	B[0]    = K[2];

	if(CRW[1])//the 3-wave is a rarefaction wave
	{
	  CSTAR2= CR * pow(P_STAR/WR[2],G1);//CSTAR2 is the sonic speed of the right star region
	  RHO2  = WR[0] * pow(P_STAR/WR[2],1.0/GAMMA);//RHO2 is the density of the right star region

	  K[2] = RAREFACTION_RIGHT(GAMMA, CSTAR2, WR, CR, SLOPER);
	  A[1][0] = 1.0;
	  A[1][1] = -1.0/(RHO2*CSTAR2);
	  B[1]    = K[2];
	}
	else//the 3-wave is a shock
	{
	//WW is the state of [rho,u,p] between the contact and the 3-wave
	  WW[0] = WR[0] *(MID[2]+MU2*WR[2])/(WR[2]+MU2* MID[2]);
	  WW[1] = MID[1];
	  WW[2] = MID[2];
	  CSTAR2= sqrt((GAMMA* MID[2])/WW[0]);//CSTAR2 is the sonic speed of the right star region
             
	  SHOCK_RIGHT(K, GAMMA, WW, CSTAR2, WR, CR, SLOPER);
	  A[1][0] = K[0];
	  A[1][1] = K[1];
	  B[1]    = K[2];
	}
	detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	K[0] = (A[1][1]*B[0] - A[0][1]*B[1]) / detA;
	K[1] = (A[0][0]*B[1] - A[1][0]*B[0]) / detA;

	PT[1] = K[0] + MID[1]/(MID[0]*CSTAR*CSTAR)* K[1];
	PT[2] = K[1] + MID[0] *MID[1] * K[0];
	PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WL, CL, SLOPEL);
      }
    }//end of the case of the left rarefaction wave
    else
    {
      QL = P_STAR/WL[2];
      SL  = WL[1] - CL*sqrt(G2*QL + G1);
      if(fabs(SL) < eps)
      {
	printf("stationary shock\n");
	getchar();
      }

      if(0.0<SL)//trivial case
      {
	MID[0] = WL[0];
	MID[1] = WL[1];
	MID[2] = WL[2];
	COEFFICIENT_MATRIX(A, GAMMA, MID);
	MATMUL(PT, A, SLOPEL);
	PT[0] = -PT[0];
	PT[1] = -PT[1];
	PT[2] = -PT[2];
      }
      else//non-trivial case
      {
	MID[2] = P_STAR;
	MID[1] = U_STAR;
	MID[0] = WL[0]*(QL + G6)/(QL * G6 + 1.0);
	CSTAR=sqrt(GAMMA * MID[2]/ MID[0]);

	SHOCK_LEFT(K, GAMMA, MID, CSTAR, WL, CL, SLOPEL);
	A[0][0] = K[0];
	A[0][1] = K[1];
	B[0]    = K[2];

        if(CRW[1])//the 3-wave is a rarefaction wave
	{
	  CSTAR2= CR * pow(P_STAR/WR[2],G1);//CSTAR2 is the sonic speed of the right star region
	  RHO2  = WR[0] * pow(P_STAR/WR[2],1.0/GAMMA);//RHO2 is the density of the right star region

	  K[2] = RAREFACTION_RIGHT(GAMMA, CSTAR2, WR, CR, SLOPER);
	  A[1][0] = 1.0;
	  A[1][1] = -1.0/(RHO2*CSTAR2);
	  B[1]    = K[2];
	}
	else//the 3-wave is a shock
	{
	  //WW is the state of [rho,u,p] between the contact and the 3-wave
	  WW[0] = WR[0] *(MID[2]+MU2*WR[2])/(WR[2]+MU2* MID[2]);
	  WW[1] = MID[1];
	  WW[2] = MID[2];
	  CSTAR2= sqrt((GAMMA* MID[2])/WW[0]);//CSTAR2 is the sonic speed of the right star region
             
	  SHOCK_RIGHT(K, GAMMA, WW, CSTAR2, WR, CR, SLOPER);
	  A[1][0] = K[0];
	  A[1][1] = K[1];
	  B[1]    = K[2];
	}
	detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	K[0] = (A[1][1]*B[0] - A[0][1]*B[1]) / detA;
	K[1] = (A[0][0]*B[1] - A[1][0]*B[0]) / detA;
      
	PT[1] = K[0] + MID[1]/(MID[0]*CSTAR*CSTAR)* K[1];
	PT[2] = K[1] + MID[0] *MID[1] * K[0];
	PT[0] = SHOCK_DENSITY_L(GAMMA, K[1], K[0], MID, CSTAR, WL, CL, SLOPEL);
      }
    }//end of the case of left shock
    UL[0] = MID[0];
    UL[1] = MID[1];
    UL[2] = MID[2];
    UR[0] = MID[0];
    UR[1] = MID[1];
    UR[2] = MID[2];
    DL[0] = PT[0];
    DL[1] = PT[1];
    DL[2] = PT[2];
    DR[0] = PT[0];
    DR[1] = PT[1];
    DR[2] = PT[2];
  }//end of the case where the t-axe lie to the left of the contact
  else if(fabs(U_STAR)<eps)
  {
    UL[1] = U_STAR;
    UL[2] = P_STAR;
    UR[1] = U_STAR;
    UR[2] = P_STAR;

    if(CRW[0])//the 1-wave is a rarefaction wave
    {
      UL[0] = WL[0] * pow(UL[2]/WL[2] , 1.0/GAMMA);
      CSTAR1 = sqrt(GAMMA*UL[2]/UL[0]);

      K[2] = RAREFACTION_LEFT(GAMMA, CSTAR1, WL, CL, SLOPEL);
      A[0][0] = 1.0;
      A[0][1] = 1.0/(UL[0]*CSTAR1);
      B[0]    = K[2];
    }
    else//the 1-wave is a shock
    {
      UL[0] = WL[0]*(UL[2]+MU2*WL[2])/(WL[2]+MU2*UL[2]);
      CSTAR1 = sqrt(GAMMA * UL[2]/ UL[0]);

      SHOCK_LEFT(K, GAMMA, UL, CSTAR1, WL, CL, SLOPEL);
      A[0][0] = K[0];
      A[0][1] = K[1];
      B[0]    = K[2];
    }
    if(CRW[1])//the 3-wave is a rarefaction wave
    {
      UR[0] = WR[0] * pow(UR[2]/WR[2] , 1.0/GAMMA);
      CSTAR2  = sqrt(GAMMA * UR[2]/ UR[0]);

      K[2] = RAREFACTION_RIGHT(GAMMA, CSTAR2, WR, CR, SLOPER);
      A[1][0] = 1.0;
      A[1][1] = -1.0/(UR[0]*CSTAR2);
      B[1]    = K[2];
    }
    else//the 3-wave is a shock
    {
      UR[0] = WR[0] * (UR[2]+MU2*WR[2]) / (WR[2]+MU2*UR[2]);
      CSTAR2= sqrt(GAMMA*UR[2]/UR[0]);

      SHOCK_RIGHT(K, GAMMA, UR, CSTAR2, WR, CR, SLOPER);
      A[1][0] = K[0];
      A[1][1] = K[1];
      B[1]    = K[2];
    }

    detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    K[0] = (A[1][1]*B[0] - A[0][1]*B[1]) / detA;
    K[1] = (A[0][0]*B[1] - A[1][0]*B[0]) / detA;

    //PT[1] = K[0] + UL[1]/(UL[0]*CSTAR1*CSTAR1)* K[1];
    //PT[1] = K[0] + UR[1]/(UR[0]*CSTAR2*CSTAR2)* K[1];
    //but, UL[1] = UR[1] = USTAR = 0
    PT[1] = K[0];
    //PT[2] = K[1] + UL[0] *UL[1] * K[0];
    //PT[2] = K[1] + UR[0] *UR[1] * K[0];
    //but, UL[1] = UR[1] = USTAR = 0
    PT[2] = K[1];

    if(CRW[0])//1-rarefaction
      DL[0] = RAREFACTION_DENSITY(GAMMA, PT[2], UL, CSTAR1, WL, CL, SLOPEL);
    else//1-shock
      DL[0] = SHOCK_DENSITY_L(GAMMA, K[1], K[0], UL, CSTAR1, WL, CL, SLOPEL);

    if(CRW[1])
      DR[0] = RAREFACTION_DENSITY(GAMMA, PT[2], UR, CSTAR2, WR, CR, SLOPER);
    else
      DR[0] = SHOCK_DENSITY_R(GAMMA, K[1], K[0], UR, CSTAR2, WR, CR, SLOPER);

    DL[1] = PT[1];
    DL[2] = PT[2];
    DR[1] = PT[1];
    DR[2] = PT[2];
  }
  else//the t-axe lies to the right side of the contact
  {
    if(CRW[1])
    {
      CSTAR= CR * pow(P_STAR/WR[2] , G1);
      HEAD_R = WR[1] + CR;
      TAIL_R = U_STAR + CSTAR;

      if(HEAD_R<0.0)//right trivial case
      {
	MID[0] = WR[0];
	MID[1] = WR[1];
	MID[2] = WR[2];
        COEFFICIENT_MATRIX(A, GAMMA, MID);
        MATMUL(PT, A, SLOPER);
	PT[0] = -PT[0];
	PT[1] = -PT[1];
	PT[2] = -PT[2];
      }
      else if(TAIL_R<0.0)//right sonic case
      {
	CSTAR  = -MU2 * RIE_RR;
	MID[1] = -CSTAR;
	MID[0] =  WR[0]*pow(CSTAR/CR , 2.0*NU);
	MID[2] =  WR[2] * pow(MID[0]/WR[0] , GAMMA);

        K[2] = RAREFACTION_RIGHT(GAMMA, CSTAR, WR, CR, SLOPER);

	PT[1] = K[2];
	PT[2] = MID[0] * MID[1] * K[2];

        PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WR, CR, SLOPER);
      }
      else//non-sonic case
      {
	MID[2] = P_STAR;
	MID[1] = U_STAR;
	MID[0] = WR[0] * pow(MID[2]/WR[2] , 1.0/GAMMA);
	CSTAR  = sqrt(GAMMA * MID[2]/ MID[0]);
	K[2] = RAREFACTION_RIGHT(GAMMA, CSTAR, WR, CR, SLOPER);
	A[1][0] = 1.0;
	A[1][1] = -1.0/(MID[0]*CSTAR);
	B[1]    = K[2];

	if(CRW[0])//left rarefaction
	{
	  CSTAR1 = CL * pow(P_STAR/WL[2] , G1);
          RHO1   = WL[0] * pow(P_STAR/WL[2] , 1.0/GAMMA);

          K[2] = RAREFACTION_LEFT(GAMMA, CSTAR1, WL, CL, SLOPEL);
	  A[0][0] =1.0;
	  A[0][1] = 1.0/(RHO1*CSTAR1);
	  B[0]    = K[2];
	}
	else//left shock
	{
	  WW[0] = WL[0] *(MID[2]+MU2*WL[2])/(WL[2]+MU2* MID[2]);
	  WW[1] = MID[1];
	  WW[2] = MID[2];
	  CSTAR1= sqrt((GAMMA* MID[2])/WW[0]);

          SHOCK_LEFT(K, GAMMA, WW, CSTAR1, WL, CL, SLOPEL);
	  A[0][0] = K[0];
	  A[0][1] = K[1];
	  B[0]    = K[2];
	}
	detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	K[0] = (A[1][1]*B[0] - A[0][1]*B[1]) / detA;
	K[1] = (A[0][0]*B[1] - A[1][0]*B[0]) / detA;

	PT[1] = K[0] + MID[1]/(MID[0]*CSTAR*CSTAR)* K[1];
	PT[2] = K[1] + MID[0] *MID[1] * K[0];
        PT[0] = RAREFACTION_DENSITY(GAMMA, PT[2], MID, CSTAR, WR, CR, SLOPER);
      }
    }
    else//right shock 
    {
      QR = P_STAR/ WR[2];
      SR  = WR[1] + CR* sqrt(G2 * QR + G1);
      if(fabs(SR)<eps)
      {
	printf("stationary shock\n");
	getchar();
      }

      if(SR<0.0)
      {
	MID[0] = WR[0];
	MID[1] = WR[1];
	MID[2] = WR[2];
	COEFFICIENT_MATRIX(A, GAMMA, MID);
	MATMUL(PT, A, SLOPER);
	PT[0] = -PT[0];
	PT[1] = -PT[1];
	PT[2] = -PT[2];
      }
      else//non-trivial case
      {
	MID[2] = P_STAR;
	MID[1] = U_STAR;
	MID[0] = WR[0] * (QR + G6)/(QR * G6 + 1.0);
	CSTAR  = sqrt(GAMMA * MID[2]/ MID[0]);

	SHOCK_RIGHT(K, GAMMA, MID, CSTAR, WR, CR, SLOPER);
	A[1][0] = K[0];
	A[1][1] = K[1];
	B[1]    = K[2];

	if(CRW[0])//left rarefaction
	{
	  CSTAR1 = CL * pow(P_STAR/WL[2] , G1);
	  RHO1   = WL[0] * pow(P_STAR/WL[2] , 1.0/GAMMA);

	  K[2] = RAREFACTION_LEFT(GAMMA, CSTAR1, WL, CL, SLOPEL);
	  A[0][0] =1.0;
	  A[0][1] = 1.0/(RHO1*CSTAR1);
	  B[0]    = K[2];
	}
	else//left shock
	{
	  WW[0] = WL[0] *(MID[2]+MU2*WL[2])/(WL[2]+MU2* MID[2]);
	  WW[1] = MID[1];
	  WW[2] = MID[2];
	  CSTAR1= sqrt((GAMMA* MID[2])/WW[0]);

	  SHOCK_LEFT(K, GAMMA, WW, CSTAR1, WL, CL, SLOPEL);
	  A[0][0] = K[0];
	  A[0][1] = K[1];
	  B[0]    = K[2];
	}
	detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	K[0] = (A[1][1]*B[0] - A[0][1]*B[1]) / detA;
	K[1] = (A[0][0]*B[1] - A[1][0]*B[0]) / detA;

	PT[1] = K[0] + MID[1]/(MID[0]*CSTAR*CSTAR)* K[1];
	PT[2] = K[1] + MID[0] *MID[1] * K[0];
	PT[0] = SHOCK_DENSITY_R(GAMMA, K[1], K[0], MID, CSTAR, WR, CR, SLOPER);
      }
    }
    UL[0] = MID[0];
    UL[1] = MID[1];
    UL[2] = MID[2];
    UR[0] = MID[0];
    UR[1] = MID[1];
    UR[2] = MID[2];
    DL[0] = PT[0];
    DL[1] = PT[1];
    DL[2] = PT[2];
    DR[0] = PT[0];
    DR[1] = PT[1];
    DR[2] = PT[2];
  }
}
