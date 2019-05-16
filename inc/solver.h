#ifndef L_BUFF
#include "file_io.h"
#endif
#ifndef bit_shift
#include "file_io_local.h"
#endif

#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */



//======================fixed GRP=================================

int GRP2_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);
//======================fixed GRP-HWENO=================================
int GRP4_HWENO5_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

int GRP4_minmod_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);


//======================fixed GRP-WENO=================================
int GRP3_WENO3_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

int GRP3_WENO5_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

int GRP4_WENO5_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

int GRP5_WENO5_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);




int RF4_WENO_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

void flux_RF(double const running_info[], int const m, double const h, double const gamma,
	     double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);
void flux_RF_dual(double const running_info[], int const m, double const h, double const gamma,
		  double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);
void flux_RF_1st(int const running_info[], int const m, double const h, double const gamma,
	     double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);

int FV_WENO_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);

int LF4_WENO_fix
(option OPT, int const m, double const h,
 double rho[], double u[], double p[], runHist *runhist,
 char *add_mkdir, char *label);





