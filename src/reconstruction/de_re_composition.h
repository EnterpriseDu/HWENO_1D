


inline static void decomposition(double H, double u, double c_square, double c,
		   double gamma1, int m, int j, int l, int shift, int flag,
		   double W1[], double W2[], double W3[],
		   double Q1[], double Q2[], double Q3[])
{
  int k;
  double Qtemp, u_square = u*u;

  if(flag)
    for(k = 0; k < l; ++k)
    {
      Qtemp = W1[j-shift+k]*0.5*u_square - W2[j-shift+k]*u + W3[j-shift+k];
      Q1[k] = Qtemp*gamma1 + c*(u*W1[j-shift+k] - W2[j-shift+k]);
      Q1[k] = 0.5*Q1[k]/c_square;
      Q2[k] = Qtemp*gamma1 - c_square*W1[j-shift+k];
      Q2[k] =-Q2[k]/c_square;
      Q3[k] = Qtemp*gamma1 - c*(u*W1[j-shift+k] - W2[j-shift+k]);
      Q3[k] = 0.5*Q3[k]/c_square;
    }
  else
    for(k = 0; k < l; ++k)
    {
      Q1[k] = W1[j-shift+k];
      Q2[k] = W2[j-shift+k];
      Q3[k] = W3[j-shift+k];
    }
}

inline static void recomposition(double H, double u, double c, int flag,
		   double QL[], double QR[], double PL[], double PR[])
{
  double u_square = u*u;
  if(flag)
  {
    PL[0] = QL[0] + QL[1] + QL[2];
    PL[1] = u*PL[0] + c*(QL[2] - QL[0]);
    PL[2] = H*(QL[0]+QL[2]) + u*c*(QL[2]-QL[0]) + 0.5*u_square*QL[1];

    PR[0] = QR[0] + QR[1] + QR[2];
    PR[1] = u*PR[0] + c*(QR[2] - QR[0]);
    PR[2] = H*(QR[0]+QR[2]) + u*c*(QR[2]-QR[0]) + 0.5*u_square*QR[1];
  }
  else
  {
    PL[0] = QL[0];
    PL[1] = QL[1];
    PL[2] = QL[2];

    PR[0] = QR[0];
    PR[1] = QR[1];
    PR[2] = QR[2];
  }
}
