#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"


void Riemann_solver_Toro(double *PP, double *U, double GAMMA, double *WL, double *WR, double eps, int N)
{
  int k = 0;
  double FL, FLD, FR, FRD, POLD, P_START, CHANGE, UDIFF, CL, CR;
  double v_L, v_R, k1, k3, temp1, temp2, temp3, delta_p;
  double mu = (GAMMA-1.0) / (2.0*GAMMA), nu = (GAMMA+1.0) / (2.0*GAMMA);

  CL=sqrt(GAMMA*WL[2]/WL[0]);
  CR=sqrt(GAMMA*WR[2]/WR[0]);

  P_START = GUESSP(GAMMA, WL, WR, eps);
  POLD  = P_START;
  UDIFF = WR[1] - WL[1];


  for(k = 0; k < N; ++k)
  {
    PREFUN(&FL, &FLD, GAMMA, POLD, WL[0], WL[2], CL);
    PREFUN(&FR, &FRD, GAMMA, POLD, WR[0], WR[2], CR);
    *PP = POLD - (FL + FR + UDIFF)/(FLD + FRD);

    if(POLD > WL[2])
    {
      delta_p = POLD - WL[2];
      temp1 = 1.0 / sqrt(1.0 + nu*delta_p/WL[2]);
      temp2 = CL / GAMMA / WL[2];
      temp3 = 0.5 * temp2 * nu / WL[2];
      k1 = temp3*delta_p*pow(temp1,3.0) - temp2*temp1;
      delta_p = POLD - WL[2];
      v_L = sqrt(1.0 + nu*delta_p/WL[2]);
      v_L = delta_p * CL / GAMMA / WL[2] / v_L;
    }
    else
    {
      temp2 = CL / GAMMA / WL[2];
      temp1 = 1.0 / pow(POLD/WL[2], nu);
      k1 = -temp1 * temp2;
      v_L = pow(POLD/WL[2], mu) - 1.0;
      v_L = 2.0 * CL * v_L / (GAMMA-1.0);
    }
    //the (p,u)-tangent slope on I3 at (v_R,POLD), i.e. [du/dp](POLD)
    if(POLD > WR[2])
    {
      delta_p = POLD - WR[2];
      temp1 = 1.0 / sqrt(1.0 + nu*delta_p/WR[2]);
      temp2 = CR / GAMMA / WR[2];
      temp3 = 0.5 * temp2 * nu / WR[2];
      k3 = temp2*temp1 - temp3*delta_p*pow(temp1,3.0);
      delta_p = POLD - WR[2];
      v_R = sqrt(1.0 + nu*delta_p/WR[2]);
      v_R = delta_p * CR / GAMMA / WR[2] / v_R;
    }
    else
    {
      temp2 = CR / GAMMA / WR[2];
      temp1 = 1.0 / pow(POLD/WR[2], nu);
      k3 = temp1 * temp2;
      v_R = pow(POLD/WR[2], mu) - 1.0;
      v_R = 2.0 * CR * v_R / (GAMMA-1.0);
    }
v_L = v_L-FL;
v_R = v_R-FR;
k1=k1+FLD;
k3=k3-FRD;
if(v_L*v_L + v_R*v_R + k1*k1 + k3*k3 > eps)
{
  printf("%d\n", k);
  getchar();
}

    CHANGE = 2.0*fabs((*PP - POLD)/(*PP + POLD));
    if(CHANGE < eps)
      break;
    if(*PP < eps)
      *PP = eps;
    POLD = *PP;
  }

  *U = 0.5*(WL[1] + WR[1] + FR - FL);
}


double GUESSP(double GAMMA, double *WL, double *WR, double eps)
{
  double G1=(GAMMA-1.0)/(2.0*GAMMA), G2=(GAMMA+1.0)/(2.0*GAMMA);
  double G3=2.0*GAMMA/(GAMMA-1.0), G4=2.0/(GAMMA-1.0);
  double G5=2.0/(GAMMA+1.0), G6=(GAMMA-1.0)/(GAMMA+1.0);
  double G7 = (GAMMA-1.0)/2.0, G8 = GAMMA - 1.0;

  double CL, CR, QUSER, CUP, PPV, PMIN, PMAX, QMAX, PM, PQ, PTL, PTR, UM, GEL, GER;

  CL=sqrt(GAMMA*WL[2]/WL[0]);
  CR=sqrt(GAMMA*WR[2]/WR[0]);

  QUSER = 2.0;
  CUP  = 0.25*(WL[0] + WR[0])*(CL + CR);
  PPV  = 0.5*(WL[2] + WR[2]) + 0.5*(WL[1] - WR[1])*CUP;
  if(PPV < eps)
    PPV = eps;
  PMIN = (WL[2] < WR[2])? WL[2] : WR[2];
  PMAX = WL[2] + WR[2] - PMIN;
  QMAX = PMAX/PMIN;

  if(QMAX < QUSER && (PMIN < PPV && PPV < PMAX))//Select PVRS Riemann solver
    PM = PPV;
  else
    if(PPV < PMIN)//Select Two-Rarefaction Riemann solver
    {
      PQ  = pow(WL[2]/WR[2] , G1);
      UM  = (PQ*WL[1]/CL + WR[1]/CR +	G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
      PTL = 1.0 + G7*(WL[1] - UM)/CL;
      PTR = 1.0 + G7*(UM - WR[1])/CR;
      PM  = 0.5*(pow(WL[2]*PTL , G3) + pow(WR[2]*PTR , G3));
    }
    else//Select Two-Shock Riemann solver with PVRS as estimate
    {
      GEL = sqrt((G5/WL[0])/(G6*WL[2] + PPV));
      GER = sqrt((G5/WR[0])/(G6*WR[2] + PPV));
      PM  = (GEL*WL[2] + GER*WR[2] - (WR[1] - WL[1]))/(GEL + GER);
      if(PM < eps)
	PM = eps;
    }

  return PM;
}



/*
 *Purpose: to evaluate the pressure functions
 *FL and FR in exact Riemann solver
 *and their first derivatives
 */
void PREFUN(double *F, double *FD, double GAMMA, double P, double DK, double PK, double CK)
{
  double G1=(GAMMA-1.0)/(2.0*GAMMA), G2=(GAMMA+1.0)/(2.0*GAMMA);
  double G3=2.0*GAMMA/(GAMMA-1.0), G4=2.0/(GAMMA-1.0);
  double G5=2.0/(GAMMA+1.0), G6=(GAMMA-1.0)/(GAMMA+1.0);
  double G7 = (GAMMA-1.0)/2.0, G8 = GAMMA - 1.0;

  double PRATIO, AK, BK, QRT;

  if(P < PK)//Rarefaction wave
  {
    PRATIO = P/PK;
    *F   = G4*CK*(pow(PRATIO , G1) - 1.0);
    *FD  = (1.0/(DK*CK))*pow(PRATIO , -G2);
  }
  else//Shock wave
  {
    AK  = G5/DK;
    BK  = G6*PK;
    QRT = sqrt(AK/(BK + P));
    *F  = (P - PK)*QRT;
    *FD = (1.0 - 0.5*(P - PK)/(BK + P))*QRT;
  }
}
