/*
 * misc.c
 *
 *  Created on: Apr 21, 2009
 *      Author: guido
 */
#include "ffoptim.h"
#include <time.h>

void Error(char *x)
{
	fprintf(stderr,"\n\n* FATAL ERROR\n");
	fprintf(stderr,"%s",x);
	fprintf(stderr,"\n\n");

	exit(1);
}

void Error2(char *x, char *y)
{
	fprintf(stderr,"\n\n* FATAL ERROR \n");
	fprintf(stderr,"%s%s",x,y);
	fprintf(stderr,"\n\n");

	exit(1);
}

double **AlloDoubleMatrix(int l, int m)
{
  double **x;
  int i,j;

  x = (double **) calloc(l,sizeof(double *));
  if (!x) Error("Cannot allocate double matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(double));
      if (!(*(x+i))) Error("Cannot allocate double matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0.;
    }

  return x;
}

int *AlloInt(int n)
{
	int *x,i;
	
	x = (int *) calloc(n,sizeof(int));
	if (!x) Error("Cannot allocate int vector");

	for (i=0;i<n;i++) x[i]=0;

	return x;
}

double *AlloDouble(int n)
{
	double *x;
	int i;
	
	x = (double *) calloc(n,sizeof(double));
	if (!x) Error("Cannot allocate double vector");

	for (i=0;i<n;i++) x[i]=0.;

	return x;
}

char **AlloString(int n, int strl)
{
    char **x;
    int i;
    
    x = (char **) calloc(n,sizeof(char *));
    if (!x) Error("Cannot allocate string");

    for (i=0;i<n;i++)
    {
        *(x+i) = calloc(strl,sizeof(char));
        if (!(*(x+i))) Error("Cannot allocate string");
        *(*(x+i)+0) = '\0';
    }
    
    return x;
}

int **AlloIntMatrix(int l, int m)
{
  int **x;
  int i,j;

  x = (int **) calloc(l,sizeof(int *));
  if (!x) Error("Cannot allocate int matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(int));
      if (!(*(x+i))) Error("Cannot allocate int matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0;
    }

  return x;
}

/*****************************************************************************
 Dumb generator of integer random numbers
 *****************************************************************************/
int irand(int r)
{
  int i;
  double q;


  q=(double) rand()/(RAND_MAX+1.0);
  q *= r;
  i= (int) q;
  return i;
}

/*****************************************************************************
 Dumb generator of double random numbers
 *****************************************************************************/
double frand()
{
  double q;

  q=(double) rand()/(RAND_MAX+1.0);
  return q;
}

/*****************************************************************************
 Initialize random number with seed (-1 to use computer time)
 *****************************************************************************/
long Randomize(int n)
{
  int i;


  if (n==-1)
  {
          n =  (time(0) - 760650080) % 3600;
          fprintf(stderr,"Random seed = %d\n",n);
          return n;
  }
  for (i=0;i<n;i++) irand(100);

  return irand(10000);
}

int IsZero(double x)
{
  if (x>-EPSILON && x<EPSILON) return 1;
  else return 0;
}


