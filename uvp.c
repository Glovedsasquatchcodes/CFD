#include "uvp.h"

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
	for (i=1; i<=imax-1; i++)
	{
		for (j=1; j<=jmax; j++)
		{
			F[i][j] = U[i][j] + dt*( ( (U[i+1][j] - 2*U[i][j] + U[i-1][j])/(dx*dx) + (U[i][j+1] - 2*U[i][j] + U[i][j-1])/(dy*dy) )/Re  - ((((U[i][j] + U[i+1][j])/2)^2 - ((U[i-1][j] + U[i][j])/2)^2 )/dx + ((|U[i][j] + U[i+1][j]|)*(U[i][j] - U[i+1][j])/4 - (|U[i-1][j] + U[i][j]|)*(U[i-1][j] - U[i][j])/4)*alpha/dx) -  ( ( (V[i][j] + V[i+1][j])*(U[i][j] + U[i][j+1])/4 - (V[i][j-1] + V[i+1][j-1])*(U[i][j-1] + U[i][j])/4 )/dy + ( (|V[i][j] + V[i+1][j]|)*(U[i][j] - U[i][j+1])/4 - (|V[i][j-1] + V[i+1][j-1]|)*(U[i][j-1] - U[i][j])/4 )*alpha/dy)  +  Gx );
		}
	}
	
	for (i=1; i<=imax; i++)
	{
		for (j=1; j<=jmax-1; j++)
		{
			G[i][j] = V[i][j] + dt*( ( (V[i+1][j] - 2*V[i][j] + V[i-1][j])/(dx*dx) + (V[i][j+1] - 2*V[i][j] + V[i][j-1])/(dy*dy) )/Re  - ((((V[i][j] + V[i][j+1])/2)^2 - ((V[i][j-1] + V[i][j])/2)^2 )/dy + ((|V[i][j] + V[i][j+1]|)*(V[i][j] - V[i][j+1])/4 - (|V[i][j-1] + V[i][j]|)*(V[i][j-1] - V[i][j])/4)*alpha/dy) -  ( ( (U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j])/4 - (U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j])/4 )/dx + ( (|U[i][j] + U[i][j+1]|)*(V[i][j] - V[i+1][j])/4 - (|U[i-1][j] + U[i-1][j+1]|)*(V[i][j-1] - V[i][j])/4 )*alpha/dx)  +  Gy );
		}
	}
	
	for (int j = 1; j <= jmax; j++)
	{
		
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
			
	}
	
	for (int i = 1; i <= imax; i++)
	{
		
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
		
	}
	
}



void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
	
	for (i=1; i<=imax; i++)
	{
		for (j=1; j<=jmax; j++)
		{
			RS[i][j] = ((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j]-1)/dy)/dt;
		}
	}
	
}



void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
	
}




void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
	for (i=1; i<=imax-1; i++)
	{
		for (j=1; j<=jmax; j++)
		{
			U[i][j] = F[i][j] - dt*(P[i+1][j] - P[i][j])/dx;
		}
	}
	
	for (i=1; i<=imax; i++)
	{
		for (j=1; j<=jmax-1; j++)
		{
			V[i][j] = G[i][j] - dt*(P[i][j+1] - P[i][j])/dy;
		}
	}
}



