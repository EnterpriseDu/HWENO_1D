double Riemann_solver_exact(double * U_star, double * P_star, double gamma, double U_l, double U_r, double P_l, double P_r, double c_l, double c_r, int * CRW, double tol, int N);

void Riemann_solver_Toro(double *PP, double *U, double GAMMA, double *WL, double *WR, double eps, int N);
double GUESSP(double GAMMA, double *WL, double *WR, double eps);
void PREFUN(double *F, double *FD, double GAMMA, double P, double DK, double PK, double CK);




void linear_GRP_solver
(double *wave_speed, double *direvative, double *source, double lambda,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   v_L, double   v_R, double   s_v_L, double   s_v_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);
void vacuum
(double *wave_speed, double *direvative, double *source, double lambda,
 double rho_L, double rho_R, double d_rho_L, double d_rho_R,
 double   u_L, double   u_R, double   d_u_L, double   d_u_R,
 double   v_L, double   v_R, double   d_v_L, double   d_v_R,
 double   p_L, double   p_R, double   d_p_L, double   d_p_R,
 double gamma, double eps);

void linear_GRP_solver_Li(double *UL, double *UR, double *DL, double *DR, 
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   v_L, double   v_R, double   s_v_L, double   s_v_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);

void    VACUUM(double *UL, double *UR, double *DL, double *DR, double GAMMA, double eps, double *WL, double *WR, double *SLOPEL, double *SLOPER);

void NONVACUUM(double *UL, double *UR, double *DL, double *DR, double GAMMA, double eps, double *WL, double *WR, double *SLOPEL, double *SLOPER);


void ACOUSTIC(double *DU, double *WL, double CL, double *WR, double CR, double *SLOPEL, double *SLOPER);

/*
 * compute d_L when there is a 1-rarefaction wave
 *
 * C0 is c_star_L, i.e. the sonic speed in the left star region
 * WL is the left state of [rho,u,p]
 * CL is the sonic speed in the left refion
 * SLOPEL is the spatial direvative of the left state of [rho,u,p]
 */
double RAREFACTION_LEFT(double GAMMA, double C0, double *WL, double CL, double *SLOPEL);


/*
 * compute d_R when there is a 3-rarefaction wave
 *
 * C0 is c_star_R, i.e. the sonic speed in the right star region
 * WR is the right state of [rho,u,p]
 * CR is the sonic speed in the right refion
 * SLOPER is the spatial direvative of the right state of [rho,u,p]
 */
double RAREFACTION_RIGHT(double GAMMA, double C0, double *WR, double CR, double *SLOPER);


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
double RAREFACTION_DENSITY(double GAMMA, double DP, double *W0, double C0, double *WW, double CC, double *SLOPE);


/*
 * compute a_L, b_L and d_L when there is a 1-shock
 *
 * W0 is the state of [rho,u,p] in the left star region
 * C0 is c_star_L, i.e. the sonic speed in the left star region
 * WL is the left state of [rho,u,p]
 * CL is the sonic speed in the left region
 * SLOPEL is the spatial direvative of the left state of [rho,u,p]
 */
void SHOCK_LEFT(double *K, double GAMMA, double *W0, double C0, double *WL, double CL, double *SLOPEL);



/*
 * compute a_R, b_R and d_R when there is a 1-shock
 *
 * W0 is the state of [rho,u,p] in the right star region
 * C0 is c_star_R, i.e. the sonic speed in the right star region
 * WR is the right state of [rho,u,p]
 * CR is the sonic speed in the right region
 * SLOPER is the spatial direvative of the right state of [rho,u,p]
 */
void SHOCK_RIGHT(double *K, double GAMMA, double *W0, double C0, double *WR, double CR, double *SLOPER);



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
double SHOCK_DENSITY_L(double GAMMA, double DP, double DU, double *W0, double C0, double *WL, double CL, double *SLOPEL);


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
double SHOCK_DENSITY_R(double GAMMA, double DP, double DU, double *W0, double C0, double *WR, double CR, double *SLOPER);

void COEFFICIENT_MATRIX(double A[3][3], double GAMMA, double *WWW);

void MATMUL(double *PT, double A[3][3], double *B);
