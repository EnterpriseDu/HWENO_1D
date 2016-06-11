//#define t(j,m) ((j < m-1) ? j : m-1)
//#define b(j,m) ((j > 0) ? j : 0)
#define t(j,m) (j+m)%m
#define b(j,m) (j+m)%m


#define N_RUNNING 9


void GRP_minmod
(double const running_info[], int const m, double const h, double const alp2,
 double const rho[], double const u[], double const p[], double const rhoI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[]);

void GRP_minmod0
(double const running_info[], int const m, double const h, double const alp2,
 double const rho[], double const u[], double const p[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[]);




void WENO_5_noD
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[]);

void WENO_5
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[], double const rhoI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[]);

void WENO_50
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[]);

void WENO_5_LF
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[]);


void local_WENO_5_inter_d(double h, double Q[8], double DQQ[2]);

void local_WENO_5_interleft(double h, double Q[8]);

void local_WENO_5_interright(double h, double Q[8]);

void local_WENO_5_inter(double h, double Q[8]);


void local_WENO_5_inter_d_Z(double h, double Q[8], double DQQ[2]);

void local_WENO_5_interleft_Z(double h, double Q[8]);

void local_WENO_5_interright_Z(double h, double Q[8]);

void local_WENO_5_inter_Z(double h, double Q[8]);





void HWENO_5
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double const rhoI[], double const momI[], double const eneI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[]);


void HWENO_5_limited
(double const running_info[], int const m, double const h, double const eps, double const alp2, double const gamma,
 double const rho[], double const mom[], double const ene[],
 double const rhoI[], double const momI[], double const eneI[], double const uI[], double const pI[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[]);

void local_HWENO_5_inter_d(double h, double Q[6], double DQ[6]);

void local_HWENO_5_inter(double h, double Q[6], double DQ[4]);

void local_HWENO_5_inter_d_Z(double h, double Q[6], double DQ[6]);

void local_HWENO_5_inter_Z(double h, double Q[6], double DQ[4]);




void THINC_primitive_0
(double const running_info[], int const m, double const h, double const thickness,
 double const rho[], double const u[], double const p[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[]);
void THINC_conservative_0
(double const running_info[], int const m, double const h, double const thickness, double const gamma1,
 double const rho[], double const u[], double const p[],
 double rho_L[], double rho_R[], double u_L[], double u_R[], double p_L[], double p_R[],
 double D_rho_L[], double D_rho_R[], double D_u_L[], double D_u_R[], double D_p_L[], double D_p_R[], int trouble[]);

void THINC_local(double result[], double const u_min, double const u_jump, double const u_bar, double const thickness, double const h);
