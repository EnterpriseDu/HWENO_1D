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

#include "file_io_local.h"
#include "solver.h"
#include "reconstruction.h"




double * RHO0 = NULL;
double * U0 = NULL;
double * V0 = NULL;
double * P0 = NULL;
double * X0 = NULL;
double * Y0 = NULL;
int L_STR = 1000;
int const bit_shift = 8;



int main(int argc, char *argv[])
{
  int len_prob = 0, i = 0, l = 0, j = 0, k = 0;
  // len_prob is the number of characters in argv[1] EXCLUDING '\0'
  while(argv[1][len_prob] != '\0')
    ++len_prob;
  char addRHO[100];
  char addU[100];
  char addV[100];
  char addP[100];
  char addX[100];
  char addY[100];
  char addCONF[100];
  char addOPT[100];
  char extend[5] = ".txt\0";
  char add[100] = "../DATA/";
  for(i = 0; i < len_prob; ++i)
    add[8+i] = argv[1][i];
  add[8+len_prob] = '/';
  for(i = 0; i < len_prob+9; ++i)
  {
    addRHO[i] = add[i];
    addU[i] = add[i];
    addV[i] = add[i];
    addP[i] = add[i];
    addX[i] = add[i];
    addY[i] = add[i];
    addCONF[i] = add[i];
    addOPT[i] = add[i];
  }
  addRHO[len_prob+9] = 'R';
  addRHO[len_prob+10] = 'H';
  addRHO[len_prob+11] = 'O';
  addU[len_prob+9] = 'U';
  addV[len_prob+9] = 'V';
  addP[len_prob+9] = 'P';
  addX[len_prob+9] = 'X';
  addY[len_prob+9] = 'Y';
  addCONF[len_prob+9] = 'C';
  addCONF[len_prob+10] = 'O';
  addCONF[len_prob+11] = 'N';
  addCONF[len_prob+12] = 'F';
  addOPT[len_prob+9] = 'O';
  addOPT[len_prob+10] = 'P';
  addOPT[len_prob+11] = 'T';
  for(i = 0; i < len_prob; ++i)
  {
    addRHO[len_prob+12+i] = argv[1][i];
    addU[len_prob+10+i] = argv[1][i];
    addV[len_prob+10+i] = argv[1][i];
    addP[len_prob+10+i] = argv[1][i];
    addX[len_prob+10+i] = argv[1][i];
    addY[len_prob+10+i] = argv[1][i];
    addCONF[len_prob+13+i] = argv[1][i];
    addOPT[len_prob+12+i] = argv[1][i];
  }
  for(i = 0; i < 5; ++i)
  {
    addRHO[len_prob+12+len_prob+i] = extend[i];
    addU[len_prob+10+len_prob+i] = extend[i];
    addV[len_prob+10+len_prob+i] = extend[i];
    addP[len_prob+10+len_prob+i] = extend[i];
    addX[len_prob+10+len_prob+i] = extend[i];
    addY[len_prob+10+len_prob+i] = extend[i];
    addCONF[len_prob+13+len_prob+i] = extend[i];
    addOPT[len_prob+12+len_prob+i] = extend[i];
  }
  i = 0;

  char err_msg[L_STR];
  int read_state, err_code = 0;
  double CONFIG[N_CONF];


  read_state = configurate(CONFIG, addCONF, err_msg);
  if(read_state)
  {
    printf("%s", err_msg);
    err_code = read_state>>bit_shift;
    printf("Error code: %02x %02x\n\n", err_code, read_state-(err_code<<bit_shift));
    return read_state;
  }
  display_config(CONFIG, argv[1]);


  printf("\n");
  return 0;
}
