#include<math.h>



void GRP_minmod
(double const running_info[], int const m, double const h, double const alp2,
 double const rho[], double const u[], double const p[], double const rhoI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int i, j;
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];

  double SL, SR, Stmp;
  double Drho[m], Du[m], Dp[m];
  double P1[m+2], P2[m+2], P3[m+2];
  for(j = 0; j < m; ++j)
  {
    P1[j+1] = rho[j];
    P2[j+1] =   u[j];
    P3[j+1] =   p[j];
  }

  if(bod < 0)
  {
    P1[0]   = rho[0];  P2[0]   =-u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] =-u[m-1];P3[m+1] = p[m-1];
  }
  else if(bod > 0)
  {
    P1[0]   = rho[0];  P2[0]   = u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] = u[m-1];P3[m+1] = p[m-1];
  }
  else
  {
    P1[0]   = rho[m-1];P2[0]   = u[m-1];P3[0]   = p[m-1];
    P1[m+1] = rho[0];  P2[m+1] = u[0];  P3[m+1] = p[0];
  }

  for(j = 0; j < m; ++j)
  {
    Drho[j] = (rhoI[j+1] - rhoI[j]);
      Du[j] = (  uI[j+1] -   uI[j]);
      Dp[j] = (  pI[j+1] -   pI[j]);
  }
  for(j = 1; j < m+1; ++j)
  {
    SR = alp2*(P1[j+1] - P1[j]);
    SL = alp2*(P1[j] - P1[j-1]);
    if(SR*SL < 0.0)
      Drho[j-1] = 0.0;
    else if(SR*Drho[j-1] < 0.0)
      Drho[j-1] = 0.0;
    else
    {
      Stmp = ( (fabs(SL) < fabs(SR)) ? SL : SR );
      if(fabs(Stmp) < fabs(Drho[j-1]))
        Drho[j-1] = Stmp;
    }

    SR = alp2*(P2[j+1] - P2[j]);
    SL = alp2*(P2[j] - P2[j-1]);
    if(SR*SL < 0.0)
      Du[j-1] = 0.0;
    else if(SR*Du[j-1] < 0.0)
      Du[j-1] = 0.0;
    else
    {
      Stmp = ( (fabs(SL) < fabs(SR)) ? SL : SR );
      if(fabs(Stmp) < fabs(Du[j-1]))
        Du[j-1] = Stmp;
    }

    SR = alp2*(P3[j+1] - P3[j]);
    SL = alp2*(P3[j] - P3[j-1]);
    if(SR*SL < 0.0)
      Dp[j-1] = 0.0;
    else if(SR*Dp[j-1] < 0.0)
      Dp[j-1] = 0.0;
    else
    {
      Stmp = ( (fabs(SL) < fabs(SR)) ? SL : SR );
      if(fabs(Stmp) < fabs(Dp[j-1]))
        Dp[j-1] = Stmp;
    }
    //Dp[j-1] = 1.4*(P3[j]/P1[j])*Drho[j-1];
  }


//==========LIMIT VALUE===========
  for(j = 0; j < m; ++j)
  {
    rho_L[j+1] = rho[j] + 0.5*Drho[j];
    rho_R[j]   = rho[j] - 0.5*Drho[j];
      u_L[j+1] =   u[j] + 0.5*  Du[j];
      u_R[j]   =   u[j] - 0.5*  Du[j];
      p_L[j+1] =   p[j] + 0.5*  Dp[j];
      p_R[j]   =   p[j] - 0.5*  Dp[j];
    D_rho_L[j+1] = Drho[j]/h;
    D_rho_R[j]   = D_rho_L[j+1];
      D_u_L[j+1] =   Du[j]/h;
      D_u_R[j]   =   D_u_L[j+1];
      D_p_L[j+1] =   Dp[j]/h;
      D_p_R[j]   =   D_p_L[j+1];
  }

  if(bod < 0)
  {
  u_L[m] = 0.0;
  u_R[m-1] = 2.0*u[m-1];
  D_u_L[m] = -u_R[m-1]/h;
  D_u_R[m-1] = D_u_L[m];
  u_R[0] = 0.0;
  u_L[1] = 2.0*u[0];
  D_u_R[0] = u_L[1] / h;
  D_u_L[1] = D_u_R[0];

  rho_R[m] = rho_L[m];
    u_R[m] =  -u_L[m];
    p_R[m] =   p_L[m];
  rho_L[0] = rho_R[0];
    u_L[0] =  -u_R[0];
    p_L[0] =   p_R[0];

  D_rho_R[m] = -D_rho_L[m];
    D_u_R[m] =    D_u_L[m];
    D_p_R[m] =   -D_p_L[m];
  D_rho_L[0] = -D_rho_R[0];
    D_u_L[0] =    D_u_R[0];
    D_p_L[0] =   -D_p_R[0];
  }
  else if(bod > 0)
  {
  rho_R[m] = rho_R[m-1];
    u_R[m] =   u_R[m-1];
    p_R[m] =   p_R[m-1];
  rho_L[0] = rho_L[1];
    u_L[0] =   u_L[1];
    p_L[0] =   p_L[1];

  D_rho_R[m] = D_rho_R[m-1];
    D_u_R[m] =   D_u_R[m-1];
    D_p_R[m] =   D_p_R[m-1];
  D_rho_L[0] = D_rho_L[1];
    D_u_L[0] =   D_u_L[1];
    D_p_L[0] =   D_p_L[1];
  }
  else
  {
  rho_R[m] = rho_R[0];
    u_R[m] =   u_R[0];
    p_R[m] =   p_R[0];
  rho_L[0] = rho_L[m];
    u_L[0] =   u_L[m];
    p_L[0] =   p_L[m];

  D_rho_R[m] = D_rho_R[0];
    D_u_R[m] =   D_u_R[0];
    D_p_R[m] =   D_p_R[0];
  D_rho_L[0] = D_rho_L[m];
    D_u_L[0] =   D_u_L[m];
    D_p_L[0] =   D_p_L[m];
  }

/*   for(j = 0; j < m+1; ++j) */
/*   { */
/*     p_L[j] = pow(rho_L[j], 1.4); */
/*     p_R[j] = pow(rho_R[j], 1.4); */
/*     D_p_L[j] = 1.4*pow(rho_L[j], 0.4)*D_rho_L[j]; */
/*     D_p_R[j] = 1.4*pow(rho_R[j], 0.4)*D_rho_R[j]; */
/*   } */
}

void GRP_minmod0
(double const running_info[], int const m, double const h, double const alp2,
 double const rho[], double const u[], double const p[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[])
{
  int i, j;
  int const    K         = (int)running_info[0];
  double const time      =      running_info[1];
  int const    half      = (int)running_info[2];
  int const    bod       = (int)running_info[3];
  int const    WENOD     = (int)running_info[4];
  int const    decomp    = (int)running_info[5];
  int const    limiter   = (int)running_info[6];
  double const threshold =      running_info[7];

  double SL, SR, Stmp;
  double Drho[m], Du[m], Dp[m];
  double P1[m+2], P2[m+2], P3[m+2];
  for(j = 0; j < m; ++j)
  {
    P1[j+1] = rho[j];
    P2[j+1] =   u[j];
    P3[j+1] =   p[j];
  }

  if(bod < 0)
  {
    P1[0]   = rho[0];  P2[0]   =-u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] =-u[m-1];P3[m+1] = p[m-1];
  }
  else if(bod > 0)
  {
    P1[0]   = rho[0];  P2[0]   = u[0];  P3[0]   = p[0];
    P1[m+1] = rho[m-1];P2[m+1] = u[m-1];P3[m+1] = p[m-1];
  }
  else
  {
    P1[0]   = rho[m-1];P2[0]   = u[m-1];P3[0]   = p[m-1];
    P1[m+1] = rho[0];  P2[m+1] = u[0];  P3[m+1] = p[0];
  }


  for(j = 1; j < m+1; ++j)
  {
    SR = alp2*(P1[j+1] - P1[j]);
    SL = alp2*(P1[j] - P1[j-1]);
    if(SR*SL < 0.0)
      Drho[j-1] = 0.0;
    else
      Drho[j-1] = ( (fabs(SL) < fabs(SR)) ? SL : SR );

    SR = alp2*(P2[j+1] - P2[j]);
    SL = alp2*(P2[j] - P2[j-1]);
    if(SR*SL < 0.0)
      Du[j-1] = 0.0;
    else
      Du[j-1] = ( (fabs(SL) < fabs(SR)) ? SL : SR );

    SR = alp2*(P3[j+1] - P3[j]);
    SL = alp2*(P3[j] - P3[j-1]);
    if(SR*SL < 0.0)
      Dp[j-1] = 0.0;
    else
      Dp[j-1] = ( (fabs(SL) < fabs(SR)) ? SL : SR );
    //Dp[j-1] = 1.4*(P3[j]/P1[j])*Drho[j-1];
  }


//==========LIMIT VALUE===========
  for(j = 0; j < m; ++j)
  {
    rho_L[j+1] = rho[j] + 0.5*Drho[j];
    rho_R[j]   = rho[j] - 0.5*Drho[j];
      u_L[j+1] =   u[j] + 0.5*  Du[j];
      u_R[j]   =   u[j] - 0.5*  Du[j];
      p_L[j+1] =   p[j] + 0.5*  Dp[j];
      p_R[j]   =   p[j] - 0.5*  Dp[j];
    D_rho_L[j+1] = Drho[j]/h;
    D_rho_R[j]   = D_rho_L[j+1];
      D_u_L[j+1] =   Du[j]/h;
      D_u_R[j]   =   D_u_L[j+1];
      D_p_L[j+1] =   Dp[j]/h;
      D_p_R[j]   =   D_p_L[j+1];
  }

  if(bod < 0)
  {
  u_L[m] = 0.0;
  u_R[m-1] = 2.0*u[m-1];
  D_u_L[m] = -u_R[m-1]/h;
  D_u_R[m-1] = D_u_L[m];
  u_R[0] = 0.0;
  u_L[1] = 2.0*u[0];
  D_u_R[0] = u_L[1] / h;
  D_u_L[1] = D_u_R[0];

  rho_R[m] = rho_L[m];
    u_R[m] =  -u_L[m];
    p_R[m] =   p_L[m];
  rho_L[0] = rho_R[0];
    u_L[0] =  -u_R[0];
    p_L[0] =   p_R[0];

  D_rho_R[m] = -D_rho_L[m];
    D_u_R[m] =    D_u_L[m];
    D_p_R[m] =   -D_p_L[m];
  D_rho_L[0] = -D_rho_R[0];
    D_u_L[0] =    D_u_R[0];
    D_p_L[0] =   -D_p_R[0];
  }
  else if(bod > 0)
  {
  rho_R[m] = rho_R[m-1];
    u_R[m] =   u_R[m-1];
    p_R[m] =   p_R[m-1];
  rho_L[0] = rho_L[1];
    u_L[0] =   u_L[1];
    p_L[0] =   p_L[1];

  D_rho_R[m] = D_rho_R[m-1];
    D_u_R[m] =   D_u_R[m-1];
    D_p_R[m] =   D_p_R[m-1];
  D_rho_L[0] = D_rho_L[1];
    D_u_L[0] =   D_u_L[1];
    D_p_L[0] =   D_p_L[1];
  }
  else
  {
  rho_R[m] = rho_R[0];
    u_R[m] =   u_R[0];
    p_R[m] =   p_R[0];
  rho_L[0] = rho_L[m];
    u_L[0] =   u_L[m];
    p_L[0] =   p_L[m];

  D_rho_R[m] = D_rho_R[0];
    D_u_R[m] =   D_u_R[0];
    D_p_R[m] =   D_p_R[0];
  D_rho_L[0] = D_rho_L[m];
    D_u_L[0] =   D_u_L[m];
    D_p_L[0] =   D_p_L[m];
  }

  /* for(j = 0; j < m+1; ++j) */
  /* { */
  /*   p_L[j] = pow(rho_L[j], 1.4); */
  /*   p_R[j] = pow(rho_R[j], 1.4); */
  /*   D_p_L[j] = 1.4*pow(rho_L[j], 0.4)*D_rho_L[j]; */
  /*   D_p_R[j] = 1.4*pow(rho_R[j], 0.4)*D_rho_R[j]; */
  /* } */
}
