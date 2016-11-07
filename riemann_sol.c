/* 
 * This is a implementation of GRP scheme for 1-D
 * scalar conservation law:
 *                    u_t + f(u)_x = 0.
 *
 * 
 * The protential of reforming this implementation to solve
 * systems of equations of consevation laws is quite limited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "file_io.h"
#ifndef L_STR
#include "file_io_local.h"
#endif
#include "solver.h"
#include "reconstruction.h"






int main(int argc, char *argv[])
{
  int i = 0, l = 0, j = 0, k = 0, it;
  char err_msg[L_STR];
  int read_state, err_code = 0;
  double CONFIG[N_CONF];
  int already_read[N_CONF];
  char ITEM[N_CONF][L_STR];


/* CONFIG[0]  is the constant of the perfect gas
 * CONFIG[1]  is the CFL number
 * CONFIG[2]  is the largest value can be seen as zero
 * CONFIG[3]  is the first limiter of the slope
 * CONFIG[4]  is the second limiter of the slope
 * CONFIG[5]  is the first parameter of the monitor function
 * CONFIG[6]  is the second parameter of the monitor function
 * CONFIG[7]  is the modifier of the mesh redistribution
 * CONFIG[8]  is the tolerance of the mesh redistribution
 * CONFIG[9]  is the parameter of identifying the discontinuities
 * CONFIG[10] is the parameter of the thickness in the THINC reconstruction
 * CONFIG[11] is the maximal step to compute.
 * CONFIG[12] is the time to stop the computation
 * CONFIG[13] is the switch of whether use an adaptive mesh
 * CONFIG[14] denote the kind of boundary condition
 * CONFIG[15] indicates whether the initial data are the primitive
 *            variables [1], or the conservative ones [0]
 * CONFIG[16] indicates whether we use the smooth derivatives [0],
 *            or the WENO-type ones in the reconstruction
 * CONFIG[17] is the switch of whether use the limiter in the reconstruction
 * CONFIG[18] is the switch of whether use the characteristic decomposition
 * CONFIG[19] is the scaling
 */
/*
  read_state = configurate(CONFIG, ITEM, already_read, err_msg, "CONF\0", argv[1]);
  if(read_state)
  {
    printf("%s", err_msg);
    err_code = read_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, read_state-(err_code<<bit_shift));
    return read_state;
  }
  display_config(CONFIG, argv[1]);


  int nInitValue = 4;
  char addInitValue[nInitValue][L_STR];
  strcpy(addInitValue[0], "RHO\0");
  strcpy(addInitValue[1], "U\0\0");
  strcpy(addInitValue[2], "P\0");
  strcpy(addInitValue[3], "X\0");
  int sizeInitValue[nInitValue];
  realArray InitValue[nInitValue];
  for(it = 0; it < nInitValue; ++it)
    InitValue[it].n_box = 0;


  read_state = initialize(nInitValue, InitValue, sizeInitValue, err_msg, addInitValue, argv[1]);
  if(read_state)
  {
    printf("%s", err_msg);
    err_code = read_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, read_state-(err_code<<bit_shift));
    return read_state;
  }
  int m = sizeInitValue[0];

  double rhoL, uL, pL, rhoR, uR, pR, h, u1, u2;
  find_cargo_realArray(&rhoL, 0, &InitValue[0]);
  find_cargo_realArray(&uL, 0, &InitValue[1]);
  find_cargo_realArray(&pL, 0, &InitValue[2]);
  find_cargo_realArray(&rhoR, m-1, &InitValue[0]);
  find_cargo_realArray(&uR, m-1, &InitValue[1]);
  find_cargo_realArray(&pR, m-1, &InitValue[2]);
  find_cargo_realArray(&h, 0, &InitValue[3]);
  for(j = 0; j < m-1;++j)
  {
    find_cargo_realArray(&u1, j, &InitValue[0]);
    find_cargo_realArray(&u2, j+1, &InitValue[0]);
    if(fabs(u1-u2) > eps)
      break;
  }
  for(k = 0; k < nInitValue; ++k)
    delete_realArray(InitValue +k);
  double x0; x0 = (double)j * h;
  if(j == m-2)
  {
    printf("Wrong data!\n");
    return 0xFF;
  }
*/
  double gamma = 1.4;//CONFIG[0]
  double rhoL = 3000.0, uL = 0.0, pL = 3000.0;
  double rhoR = 9.850612054411153, uR = 4.031118586431574, pR = 1.0;
  double cL = sqrt(gamma * pL / rhoL);
  double cR = sqrt(gamma * pR / rhoR);
  double u_star, p_star;
  int CRW[2];
  Riemann_solver_exact(&u_star, &p_star, gamma, uL, uR, pL, pR, cL, cR, CRW, 1e-14, 1000);
  if(fabs(u_star) < 1e-9)
    u_star = 0.0;
  double c_star_L, c_star_R, rho_star_L, rho_star_R, speed_L, speed_R, zeta;
  if(CRW[0]){
    rho_star_L = rhoL*pow(p_star/pL, 1.0/gamma);
      c_star_L =   cL*pow(p_star/pL, 0.5*(gamma-1.0)/gamma);
       speed_L =   uL - cL;
  }  else{
    rho_star_L = rhoL*(p_star+zeta*pL)/(pL+zeta*p_star);
      c_star_L = sqrt(gamma * p_star / rho_star_L);
       speed_L = uL - cL*sqrt(0.5*((gamma+1.0)*(p_star/pL) + (gamma-1.0))/gamma);
  }
  if(CRW[1]){
    rho_star_R = rhoR*pow(p_star/pR,1.0/gamma);
      c_star_R =   cR*pow(p_star/pR, 0.5*(gamma-1.0)/gamma);
       speed_R =   uR + cR;
  }  else{
    rho_star_R = rhoR*(p_star+zeta*pR)/(pR+zeta*p_star);
      c_star_R = sqrt(gamma * p_star / rho_star_R);
       speed_R = uR + cR*sqrt(0.5*((gamma+1.0)*(p_star/pR) + (gamma-1.0))/gamma);
  }
  printf("%g \t | \t %g \t | \t %g \t | \t %g \n", rhoL, rho_star_L, rho_star_R, rhoR);
  printf("%g \t | \t %g \t | \t %g \t | \t %g \n", uL, u_star, u_star, uR);
  printf("%g \t | \t %g \t | \t %g \t | \t %g \n", pL, p_star, p_star, pR);
  printf("%g \t | \t %g \t | \t %g \t | \t %g \n", cL, c_star_L, c_star_R, cR);
  printf("%d, %d\n", CRW[0], CRW[1]);


  /*
  char version[L_STR] = "\0", scheme[L_STR] = "exact\0", err_msg[L_STR];
  char add_mkdir[L_STR+L_STR];
  strcpy(add_mkdir, "../SOLUTION/\0");
  CONFIG[13] = 0;
  CONFIG[16] = 0;
  CONFIG[18] = 1;
  CONFIG[17] = 0;
  state = make_directory(add_mkdir, err_msg, argv[2], scheme, version, 0, 0, CONFIG);
  if(state)
  {
    printf("%s", err_msg);
    exit(state);
  }


  int output_flag[nInitValue], output_idx[nInitValue], output_state;
  output_flag[0] = 1; //rho
  output_flag[1] = 1; //u
  output_flag[2] = 1; //p
  output_flag[3] = 1; //x
  for(it = 0; it < nInitValue; ++it)
    output_idx[it] = 0;
  
  output_state = DATA_OUTPUT(err_msg, nInitValue, addInitValue, sizeInitValue, DATA_OUT, output_idx, output_flag, add_mkdir);
  if(output_state)
  {
    printf("%s", err_msg);
    err_code = output_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, output_state-(err_code<<bit_shift));
    delete_runHist(&runhist);
    return output_state;
  }

  

  delete_runHist(&runhist);
  */
  printf("\n");
  return 0;
}
