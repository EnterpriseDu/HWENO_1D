#include <math.h>
#include <stdlib.h>


void area_calculator
(double area_norm[4], double tau, int part,
 double x1, double x2, double x3, double x4,
 double y1, double y2, double y3, double y4)
{
  area_norm[0] = 0.0;
  area_norm[1] = 0.0;
  area_norm[2] = 0.0;
  area_norm[3] = 0.0;

  int i, j;

  double h = 1.0/(double)part, a, b;
  double a1 = x4-x1, a2 = x3-x2, b1 = y4-y1, b2 = y3-y2;

  double A = (a1-a2)*(a1-a2) + (b1-b2)*(b1-b2);
  double B = a2*(a1-a2) + b2*(b1-b2);
  double C = (x1-x2)*(a1-a2) + (y1-y2)*(b1-b2);
  double D1 = a2*a2 + b2*b2 + tau*tau;
  double D2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
  double D3 = a2*(x1-x2) + b2*(y1-y2);
  double E, F, G;

  double nx, ny, nt, norm, Guass;
  for(j = 0; j < part; ++j)
    for(i = 0; i < part; ++i)
    {
      a = (0.5+(double)j)*h;
      b = (0.5+(double)i)*h;
      E = A*b*b + 2.0*B*b + D1;
      G = A*a*a + 2.0*C*a + D2;
      F = A*a*b + B*a + C*b + D3;
      Guass = sqrt(E*G - F*F);

      nx = tau*((y1-y2) + a*(b1-b2));
      ny = tau*((x2-x1) + a*(a2-a1));
      nt = a*(a1*b2 - a2*b1) + b*((x1-x2)*(b1-b2) + (y1-y2)*(a1-a2)) + b2*(x1-x2) - a2*(y1-y2);
      norm = sqrt(nx*nx + ny*ny + nt*nt);
      nx = nx/norm;
      ny = ny/norm;
      nt = nt/norm;

      area_norm[0] += Guass;
      area_norm[1] += nx*Guass;
      area_norm[2] += ny*Guass;
      area_norm[3] += nt*Guass;
    }

  area_norm[0] = area_norm[0]*h*h;
  area_norm[1] = area_norm[1]*h*h;
  area_norm[2] = area_norm[2]*h*h;
  area_norm[3] = area_norm[3]*h*h;

  area_norm[1] = area_norm[1] / area_norm[0];
  area_norm[2] = area_norm[2] / area_norm[0];
  area_norm[3] = area_norm[3] / area_norm[0];
}
