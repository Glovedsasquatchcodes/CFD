#include <helper.h>>
#include <visual.h>
#include <init.h>
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */	
int main(int argn, char** args){
					double xlength;           /* length of the domain x-dir.*/
                    double ylength;           /* length of the domain y-dir.*/
                    int  imax;                /* number of cells x-direction*/
                    int  jmax;                /* number of cells y-direction*/
                    double dx;                /* length of a cell x-dir. */
                    double dy;                /* length of a cell y-dir. */
                    double t ;              /* gravitation y-direction */
                    double t_end;             /* end time */
                	double dt;                /* time step */
                	double tau;               /* safety factor for time step*/
                	double dt_value;          /* time interval for writing visualization data in afile*/
                	int n;                    /* current time iteration step*/
                	int itermax;             /* max. number of iterations  */
                    int it;                   /* SOR iteration counter*/
                    double res;                /*residual norm of the pressure equation*/
                    double eps;               /* accuracy bound for pressure*/
                    double alpha;             /* uppwind differencing factor*/
                    double omg;               /* relaxation factor */
					
                    double Re;                /* reynolds number   */
                    double UI;                /* velocity x-direction */
                    double VI;                /* velocity y-direction */
                    double PI;                /* pressure */
                    double GX;                /* gravitation x-direction */
                    double GY;
					const char *szFileName;       /* name of the file */
                    
		    
	 read_parameters(  &szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value); 
		    
		    double **U = matrix(1,imax-1,1,jmax);
		    double **V = matrix(1,imax,1,jmax-1);
		    double **P = matrix(1,imax,1,jmax);
		 	double **RS = matrix(1,imax,1,jmax);
		 	double **F = matrix(1,imax-1,1,jmax);
		 	double **G = matrix(1,imax,1,jmax-1);
		 	
		 	t=0;
			n=0;
			
			init_uvp(UI,VI,PI,imax,jmax,&U,&V,&P);

			while t<t_end
			{
				calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,&U,&V);
				boundaryvalues(imax,jmax,&U,&V);
				calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,&U,&V,&F,&G);
				calculate_rs(dt,dx,dy,imax,jmax,&F,&G,&RS);
				
				it=0;
				while it<itermax && res>eps
				{
					sor(omg,dx,dy,imax,jmax,&P,&RS,&res);
					it=it+1;
				}
				
				calculate_uv(dt,dx,dy,imax,jmax,&U,&V,&F,&G,&P);
				
				t=t+dt;
				n=n+1;
			}
  return -1;
}
