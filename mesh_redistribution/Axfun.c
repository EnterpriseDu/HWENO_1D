void Axfun(int m, int n, double * Ap, double * OMEGA[m], double * p)
{
  int i = 0, j = 0;
  double diag, add;

  for(j = 1; j < m; ++j)
    for(i = 1; i < n; ++i)
    {
      diag = OMEGA[j][i] + OMEGA[j][i-1] + OMEGA[j-1][i] + OMEGA[j-1][i-1];
      add = 0.0;

      Ap[j*n -j -n +i] = p[j*n -j -n +i];

      if(j > 1)
        add -= 0.5*(OMEGA[j-1][i]+OMEGA[j-1][i-1]) * p[j*n -n -n -j +i +1];
      if(j < m-1)
        add -= 0.5*(OMEGA[j][i]+OMEGA[j][i-1]) * p[j*n -j +i -1];
      if(i > 1)
        add -= 0.5*(OMEGA[j-1][i-1]+OMEGA[j][i-1]) * p[j*n -n -j +i -1];
      if(i < n-1)
        add -= 0.5*(OMEGA[j-1][i]+OMEGA[j][i]) * p[j*n -n -j +i +1];

      add /= diag;
      Ap[j*n -j -n +i] += add;
    }

}
