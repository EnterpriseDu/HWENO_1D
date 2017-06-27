//#include "type_define.h"


#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */



/*
 * This function redistribute the 2-D quadrangular mesh
 * by using Winslow's variation method.
 *
 * As is well known, Winslow's variation method request
 * to solve a system of decoupled eliptic equations:
 *
 *     (omega * x_xi)_xi + (omega * x_eta)_eta = 0
 *     (omega * y_xi)_xi + (omega * y_eta)_eta = 0,
 *
 * where {xi, eta} is the logic coordinates. Each of the
 *  above equations can be discretized into a linear system
 * of equations with a five-diagnal coefficient matrix A
 * by five-point difference method.
 *
 * In this function, we use pre-conditioner conjugate
 * gradient (PCG) method, with the pre-conditioner B to
 * be a diagnal matrix whose diagnal entries is the inverse
 * of the diagnal entries of A.
 */
int redistribution
(int m, int n, double * x[m+1],   double * y[m+1], double * OMEGA[m],
 int max_it, double tol, double eps, double * number, char * charact);


void sub_redistribution
(int m, int n, double * x[m+1], double * y[m+1], double * OMEGA[m],
 int max_it, double tol, double eps, double * number, char * charact);


void Axfun(int m, int n, double * Ap, double * OMEGA[m], double * p);
void Ayfun(int m, int n, double * Ap, double * OMEGA[m], double * p);



void area_calculator
(double area_norm[4], double tau, int part,
 double x1, double x2, double x3, double x4,
 double y1, double y2, double y3, double y4);
