void Ayfun(int m, int n, double * Ap, double * OMEGA[m], double * p)
{
  int i = 0, j = 0;
  double diag, add;

  for(i = 1; i < n; ++i)
    for(j = 1; j < m; ++j)
    {
      diag = OMEGA[j][i] + OMEGA[j][i-1] + OMEGA[j-1][i] + OMEGA[j-1][i-1];
      add = 0.0;

      Ap[i*m -i -m +j] = p[i*m -i -m +j];

      if(i > 1)
        add -= 0.5*(OMEGA[j-1][i-1]+OMEGA[j][i-1]) * p[i*m -m -m -i +j +1];
      if(i < n-1)
        add -= 0.5*(OMEGA[j-1][i]+OMEGA[j][i]) * p[i*m -i +j -1];
      if(j > 1)
        add -= 0.5*(OMEGA[j-1][i]+OMEGA[j-1][i-1]) * p[i*m -i -m +j -1];
      if(j < m-1)
        add -= 0.5*(OMEGA[j][i]+OMEGA[j][i-1]) * p[i*m -i -m +j +1];

      add /= diag;
      Ap[i*m -i -m +j] += add;
    }

}
