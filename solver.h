


#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */


//======================fixed GRP=================================
int GRP2_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int ADER2_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int THINC_primitive_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int THINC_conservative_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);
//======================fixed GRPWENO=================================
int ADER4_HWENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int GRP4_HWENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int GRP4_minmod_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);



int RF4_WENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int FD_1st_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

void flux_RF(double const running_info[], int const m, double const h, double const gamma,
	     double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);

void flux_RF_dual(double const running_info[], int const m, double const h, double const gamma,
		  double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);

void flux_RF_1st(double const running_info[], int const m, double const h, double const gamma,
		 double const rho[], double const mom[], double const ene[], double F1[], double F2[], double F3[]);

int FV_WENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);

int LF4_WENO_fix
(double const CONFIG[], double const OPT[], int const m, double const h,
 double *rho[], double *u[], double *p[], runList * runhist, char *scheme);





