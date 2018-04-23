/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdio.h>
#include <math.h>
#include <uvp.h>


void calculate_dt(Re,tau,dt,dx,dy,imax,jmax,U,V){
	double U1=fabs(U[0][0]);
	double V1=fabs(V[0][0]);
	
    for(int c=0 ; c <=imax ; c++ ){
      for(int d= 0 ; d <=jmax ; d++ ){
         if ( fabs(U[c][d]) > fabs(U1) )
            U1= U[c][d];
      }
   }
    for(int c=0 ; c<=imax ; c++ ){
          for(int d=0 ; d<=jmax ; d++ ) {
             if ( fabs(U[c][d]) > fabs(U1) )
                V1= V[c][d];
          }
       }
    temp1=((Re/2)*(pow(2,dx))*(pow(2,dy))/(pow(2,dx)+pow(2,dy));
    temp2=(dx)/fabx(U);
    temp3=(dy)/fabs(V);
    temp4=fmin(temp1,temp2);
    dt=tau*fmin(temp4,temp3);
    //return dt;
}



void calculate fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G){
	double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;
	for (int j=0;j<=jmax;j++){
		F[0][j]=U[0][j];     //Boundary conditions
		F[imax][j]=U[imax][j];
	}
	
	for(int i=0;i<=imax;i++){
		G[i][0]=V[i][0];    // Boundary conditions
		G[i][jmax]=V[i][jmax]; 
	}
	
	for(int i=1; i<imax; i++){
		for (int j=1; j<jmax; j++){
			du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
			du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
			a=(U[i][j]+U[i+1][j])/2;
			b=(U[i-1][j]+U[i][j])/2;
			du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
			duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]-V[i+1][j-1])*(U[i][j-1]-U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
			F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX);
			
			dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
			dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
			c=(V[i][j]+V[i][j+1])/2;
			d=(V[i][j-1]+V[i][j])/2;
			dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
			duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]-U[i-1][j+1])*(V[i-1][j]-V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
			G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY);
			
			
		}
	}
}
