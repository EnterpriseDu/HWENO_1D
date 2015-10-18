#include <math.h>


void ACOUSTIC(double *DU, double *WL, double CL, double *WR, double CR, double *SLOPEL, double *SLOPER)
{
  double RIEL, RIER, AVER;

  RIEL = SLOPEL[1] + SLOPEL[2]/(WL[0]*CL);
  RIER = SLOPER[1] - SLOPER[2]/(WR[0]*CR);

  DU[1] = -0.5 * ( (WL[1] +CL)* RIEL + (WR[1]-CR) * RIER);
  DU[2] = -0.5 * WL[0] *CL*  ( (WL[1] +CL)* RIEL - (WR[1]-CR) * RIER) ;

  AVER = 0.5 *(WL[1] + WR[1]);

  if(AVER >0.0)
    DU[0] = (DU[2] + AVER*(SLOPEL[2]- CL*CL*SLOPEL[0])) /CL/CL;//(DU[1] + AVER*(SLOPEL[2]- CL*CL*SLOPEL[0])) /CL/CL;
  else
    DU[0] = (DU[2] + AVER*(SLOPER[2]- CR*CR*SLOPER[0])) /CR/CR;//(DU[1] + AVER*(SLOPER[2]- CR*CR*SLOPER[0])) /CR/CR
}

/*
 * compute d_L when there is a 1-rarefaction wave
 *
 * C0 is c_star_L, i.e. the sonic speed in the left star region
 * WL is the left state of [rho,u,p]
 * CL is the sonic speed in the left refion
 * SLOPEL is the spatial direvative of the left state of [rho,u,p]
 */
double RAREFACTION_LEFT(double GAMMA, double C0, double *WL, double CL, double *SLOPEL)
{
  double KK, NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, S_SLOPE_L;

  //TDS = T_L S^{\prinme}_L
  TDS = NU/WL[0]*(SLOPEL[2]-CL*CL*SLOPEL[0]);
  //S_SLOPE_L = \Psi^{\prime}_L
  S_SLOPE_L = SLOPEL[1]+NU/(WL[0]*CL)*(GAMMA*SLOPEL[2]-CL*CL*SLOPEL[0]);

  VAR1 = C0/CL;
  VAR2 = (1.0+MU2)/(1+2.0*MU2);

  //KK = d_L in paper
  KK = (VAR2* pow(VAR1,0.5/MU2) +(1.0-VAR2)*pow(VAR1,(1.0+MU2)/MU2) )*TDS - CL*pow(VAR1,0.5/MU2)*S_SLOPE_L;

  return KK;
}


/*
 * compute d_R when there is a 3-rarefaction wave
 *
 * C0 is c_star_R, i.e. the sonic speed in the right star region
 * WR is the right state of [rho,u,p]
 * CR is the sonic speed in the right refion
 * SLOPER is the spatial direvative of the right state of [rho,u,p]
 */
double RAREFACTION_RIGHT(double GAMMA, double C0, double *WR, double CR, double *SLOPER)
{
  double KK, NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, R_SLOPE_R;
 
  TDS= NU/WR[0]*(SLOPER[2]-CR*CR*SLOPER[0]);

  R_SLOPE_R= SLOPER[1]-NU/(WR[0]*CR)*(GAMMA*SLOPER[2]-CR*CR*SLOPER[0]);

  VAR1=C0/CR;
  VAR2= (1.0+MU2)/(1+2.0*MU2);
   
  KK = (VAR2* pow(VAR1,0.5/MU2) +(1.0-VAR2)*pow(VAR1,(1.0+MU2)/MU2) )*TDS +CR*pow(VAR1,0.5/MU2)*R_SLOPE_R;

  return KK;
}


/*
 * compute the temporal direvative of the density when the t-axe
 * lies between the contact and a CRW, whether it is a 1-CRW or
 * 3-CRW
 *
 * DP    is the temporal direvative of the pressure
 * W0    is the state of [rho,u,p] in the left/right star region
 *         where the t-axe lies
 * C0    is the sonic speed in the left/right star region where
 *         the t-axe lies
 * WW    is the left/right state of [rho,u,p] depends on balabala
 * CC    is the left/right sonic speed depends on balabala
 * SLOPE is the spatial direvative of WW
 */
double RAREFACTION_DENSITY(double GAMMA, double DP, double *W0, double C0, double *WW, double CC, double *SLOPE)
{
  double DD, NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR;

  TDS = NU/WW[0]*(SLOPE[2]-CC*CC*SLOPE[0]);
  VAR = pow(C0/CC , (1+MU2)/MU2);

  DD = (DP+ (GAMMA-1.0)* W0[0] *W0[1] * VAR * TDS)/C0/C0;
  return DD;
}


/*
 * compute a_L, b_L and d_L when there is a 1-shock
 *
 * W0 is the state of [rho,u,p] in the left star region
 * C0 is c_star_L, i.e. the sonic speed in the left star region
 * WL is the left state of [rho,u,p]
 * CL is the sonic speed in the left region
 * SLOPEL is the spatial direvative of the left state of [rho,u,p]
 */
void SHOCK_LEFT(double *K, double GAMMA, double *W0, double C0, double *WL, double CL, double *SLOPEL)
{
  double NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, S_SLOPE_R;

  double GES, GESL, VAR, PHIP, PHIPL, PHIDL, LP, LU, LD;

  GES = -sqrt(0.5*((GAMMA+1)*WL[2]+(GAMMA-1)*W0[2])/W0[0]);
  GESL= -sqrt(0.5*((GAMMA+1)*W0[2]+(GAMMA-1)*WL[2])/WL[0]);

  VAR = sqrt((1-MU2)/(WL[0]*(W0[2]+MU2*WL[2])));

  PHIP =  0.5* VAR * (W0[2]+(1+2.0*MU2)*WL[2])/(W0[2]+MU2*WL[2]);
  PHIPL = -0.5*VAR * ((2.0+MU2)*W0[2] + MU2*WL[2])/(W0[2]+MU2*WL[2]);
  PHIDL = -0.5*(W0[2]-WL[2])* VAR /WL[0];

  LP = -1.0/WL[0] - GESL * PHIPL;
  LU = GESL +WL[0] * (CL*CL * PHIPL + PHIDL);
  LD = -GESL * PHIDL;

 	
  K[0] = 1.0 -W0[0]* GES * PHIP;
  K[1] = -GES/(W0[0]*C0*C0)+ PHIP;
  K[2] = LP* SLOPEL[2] + LU *SLOPEL[1] + LD* SLOPEL[0];
}



/*
 * compute a_R, b_R and d_R when there is a 1-shock
 *
 * W0 is the state of [rho,u,p] in the right star region
 * C0 is c_star_R, i.e. the sonic speed in the right star region
 * WR is the right state of [rho,u,p]
 * CR is the sonic speed in the right region
 * SLOPER is the spatial direvative of the right state of [rho,u,p]
 */
void SHOCK_RIGHT(double *K, double GAMMA, double *W0, double C0, double *WR, double CR, double *SLOPER)
{
  double NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, S_SLOPE_R;

  double GES, GESR, VAR, PHIP, PHIPR, PHIDR, LP, LU, LD;


  GES  = sqrt(0.5*((GAMMA+1)*WR[2]+(GAMMA-1)*W0[2])/W0[0]);
  GESR = sqrt(0.5*((GAMMA+1)*W0[2]+(GAMMA-1)*WR[2])/WR[0]);

  VAR  = sqrt((1-MU2)/(WR[0]*(W0[2]+MU2*WR[2])));

  PHIP = 0.5* VAR * (W0[2]+(1+2.0*MU2)*WR[2])/(W0[2]+MU2*WR[2]);
  PHIPR = -0.5*VAR * ((2.0+MU2)*W0[2]+MU2*WR[2])/(W0[2]+MU2*WR[2]);
  PHIDR = -0.5*(W0[2]-WR[2])* VAR /WR[0];

  LP = -1.0/WR[0] + GESR* PHIPR;
  LU = GESR -WR[0] * (CR*CR * PHIPR + PHIDR);
  LD = GESR * PHIDR;

  K[0] = 1.0 +W0[0]* GES * PHIP;
  K[1] = -( GES/(W0[0]*C0*C0)+ PHIP);
  K[2] = LP* SLOPER[2] + LU *SLOPER[1] + LD* SLOPER[0];
}



/*
 * compute the temporal direvative of the density when the t-axe
 * lies between the contact and the 1-shock
 *
 * DP     is the TOTAL direvative of the pressure
 * DU     is the TOTAL direvative of the velocity
 * W0     is the state of [rho,u,p] in the left star region
 * C0     is the sonic speed in the left star region
 * WL     is the left state of [rho,u,p]
 * CL     is the leftt sonic speed
 * SLOPEL is the spatial direvative of WL
 */
double SHOCK_DENSITY_L(double GAMMA, double DP, double DU, double *W0, double C0, double *WL, double CL, double *SLOPEL)
{
  double DD, NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, S_SLOPE_R;

  double GES, GESL, VAR, HP, HPL, HDL, LP, LU, LD;


  GES = -sqrt(0.5*((GAMMA+1)*WL[2]+(GAMMA-1)*W0[2])/W0[0]);
  GESL = -sqrt(0.5*((GAMMA+1)*W0[2]+(GAMMA-1)*WL[2])/WL[0]);

  VAR = WL[2] + MU2* W0[2];

  HP  =  WL[0] * (1.0-MU2*MU2)*WL[2] / VAR/VAR;
  HPL = -WL[0] * (1.0-MU2*MU2)*W0[2] / VAR/VAR;
  HDL = (W0[2] +MU2* WL[2])/VAR;

  LP = HPL * GESL * SLOPEL[2];
  LD = GESL * HDL * SLOPEL[0];
  LU = -WL[0]*(HPL* CL*CL + HDL) * SLOPEL[1];

  DD = ((W0[1]+GES)/C0/C0 - HP* W0[1]) *DP +W0[0] *W0[1] * GES * HP *DU;
  DD = (DD- W0[1]*(LP+ LD +LU))/ GES;
  return DD;
}


/*
 * compute the temporal direvative of the density when the t-axe
 * lies between the contact and the 3-shock
 *
 * DP     is the TOTAL direvative of the pressure
 * DU     is the TOTAL direvative of the velocity
 * W0     is the state of [rho,u,p] in the right star region
 * C0     is the sonic speed in the right star region
 * WR     is the right state of [rho,u,p]
 * CR     is the right sonic speed
 * SLOPER is the spatial direvative of WR
 */
double SHOCK_DENSITY_R(double GAMMA, double DP, double DU, double *W0, double C0, double *WR, double CR, double *SLOPER)
{
  double DD, NU = 1.0/(GAMMA-1.0), MU2 = (GAMMA-1.0)/(GAMMA+1.0);
  double TDS, VAR1, VAR2, S_SLOPE_R;

  double GES, GESR, VAR, HP, HPR, HDR, LP, LU, LD;


  GES = sqrt(0.5*((GAMMA+1)*WR[2]+(GAMMA-1)*W0[2])/W0[0]);
  GESR = sqrt(0.5*((GAMMA+1)*W0[2]+(GAMMA-1)*WR[2])/WR[0]);

  VAR = WR[2] + MU2* W0[2];
  HP  = WR[0] * (1.0-MU2*MU2)*WR[2] / VAR/VAR;
  HPR = -WR[0] * (1.0-MU2*MU2)*W0[2] / VAR/VAR;
  HDR = (W0[2] +MU2* WR[2])/VAR;

  LP = HPR * GESR * SLOPER[2];
  LD = GESR * HDR * SLOPER[0];
  LU = -WR[0]*( HPR* CR*CR + HDR) * SLOPER[1];

  DD = ((W0[1]+GES)/C0/C0 - W0[1]*HP) *DP +W0[0]* W0[1] * GES * HP *DU;
  DD = (DD-   W0[1]*(LP+ LD +LU)   )/ GES;

  return DD;
}

void COEFFICIENT_MATRIX(double A[3][3], double GAMMA, double *WWW)
{
  double CC;
  CC=sqrt(GAMMA*WWW[2]/WWW[0]);

  A[0][0] = WWW[1]; A[0][1] = WWW[0];        A[0][2]= 0.0;
  A[1][0] = 0.0;    A[1][1] = WWW[1];        A[1][2]= 1.0/WWW[0];
  A[2][0] = 0.0;    A[2][1] = WWW[0]*CC*CC;  A[2][2]= WWW[1];
}

void MATMUL(double *PT, double A[3][3], double *B)
{
  int k;
  for(k=0; k < 3; ++k)
    PT[k] = A[k][0]*B[0] + A[k][1]*B[1] + A[k][2]*B[2];
}
