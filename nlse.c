/*********************************************************************************
* Numerical Solution for the Cubic Nonlinear Schrodinger Equation in (1+1)D	 	  *
* using symmetric split step Fourier method		                              	  *
* Coded by: Omar Ashour, Texas A&M University at Qatar, February 2015.    	      *
* ********************************************************************************/
#include "lib.h"
#define ENERGY 0 

int main(int argc, char *argv[])
{
    // Grid parameters
    double dt, A1, q, tm; int order, nx; char type;
    printf("Enter dt: ");
    scanf("%lf", &dt);
    printf("Enter tm: ");
    scanf("%lf", &tm);
    printf("Enter nx: ");
    scanf("%d", &nx);
    printf("Enter A1: "); 
    scanf("%lf", &A1);
    printf("Enter q (0<q<0.5): ");
    scanf("%lf", &q);
    printf("Enter algorithm order (2, 4, 6, 8): ");
    scanf("%d", &order);
    printf("Enter algorithm type (S for symplect, M for multiproduct): ");
    scanf(" %c", &type);

    // const double dt = atof(argv[1]);			// temporal step size
    const int nt = tm/dt;   			        // number of temporal nodes
    const double l = M_PI/sqrt(1-2*q);		    // Spatial Period
    const double dx = (l / nx);			        // spatial step size
    const double Omega = 2*sqrt(1-2*q);
    const double A0 = sqrt(1-2*A1*A1);
    
    char paramf[20];
    sprintf(paramf, "output/param.txt");
    FILE *fp = fopen(paramf, "w");
    fprintf(fp, "%f\n%d\n%f\n%f\n%f\n%d\n%c", dt, nx, tm, A1, q, 
                                              order, type);
    fclose(fp);

    // Print basic info about simulation
    printf("dt=%f\nnx=%d\ntm=%f\norder=%d\ntype=%c\nq=%f\n", 
            dt, nx, tm, order, type, q);

    // Allocate the arrays
    double *x = (double*)malloc(sizeof(double) * nx);
    double *k2 = (double*)malloc(sizeof(double) * nx);
    fftw_complex *psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    double *psi_r = (double*) malloc(sizeof(double) * nx * nt);
    double *psi_i = (double*) malloc(sizeof(double) * nx * nt);

    fftw_complex *psi1, *psi2, *psi3, *psi4;
    // Allocating extra arrays for high orders
    if (order > 2)
    {
        psi1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
        psi2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    }
    if (order > 4)
        psi3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    if (order > 6)
        psi4 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);

    // Create transform plans
    printf("Creating plans.\n");
    fftw_plan forward, backward;
    forward = fftw_plan_dft_1d(nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_1d(nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Create wave number. Flipped to make up for not shifting 0 frequency component to the middle.
    printf("Setting up initial conditions.\n");
    double dk=2*M_PI/nx/dx;
    double *k = (double*)malloc(nx*sizeof(double));
    for(int i = nx/2; i >= 0; i--) 
        k[nx/2-i]=(nx/2-i)*dk;
    for(int i = nx/2+1; i < nx; i++)
	    k[i]=(i-nx)*dk; 

    // Initial conditions
    for (int i = 0; i < nx; i++)
    {
        x[i] = (i-nx/2)*dx;
	    k2[i] = k[i]*k[i];
	    psi[i] = (A0+2*A1*cos(Omega*x[i])) + 0*I; // Don't forget 2*A1 for q=3/8
    }
   
    #if ENERGY
    double *temp1 = (double*) malloc(sizeof(double) * nx);
    double *temp2 = (double*) malloc(sizeof(double) * nx);
    
    const double VA = -1.0/2.0*(1+8*A*A-14*A*A*A*A);
    const double TA = A*A;
    double V, T;
    double s1=0, s2=0;

    for(int i = 0; i < nx; i++)
    {
        temp1[i] = psi[i]*psi[i];
        temp2[i] = temp1[i]*temp1[i];
        s1 += temp1[i];
        s2 += temp2[i];
    }
    V = -0.5*s2/s1;

    s1=0, s2=0;
    fftw_execute(forward);
    for(int i = 0; i < nx; i++)
    {
        temp1[i] = 0.5*k[i]*k[i]*psi[i]*psi[i];
        temp2[i] = psi[i]*psi[i];
        s1 += temp1[i];
        s2 += temp2[i];
    }
    T = s1/s2;
    
    char temp0[30];
    sprintf(temp0, "energy_A%0.0f.txt", -log10(A));
    fp = fopen(temp0, "w");
    fprintf(fp, "<V>+1/2 (A):%0.30e\n", VA+0.5);
    fprintf(fp, "<V>+1/2 (P):%0.30e\n", V+0.5);
    fprintf(fp, "<T> (A):%0.30e\n", TA);
    fprintf(fp, "<T> (P):%0.30e\n", T);
    fclose(fp);
    return 0;
    #endif //ENERGY
   
    // Start time evolution
    printf("Evolving time.\n");
    for (int i = 0; i < nt; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            psi_i[ind(i,j)] = cimag(psi[j]);
            psi_r[ind(i,j)] = creal(psi[j]);
        }	
 	    if (order == 2)
	        T2(psi, psi, nt, dt, k2, nx, forward, backward); 
        else if (order == 4 && type == 'M')
            T4_M(psi, psi1, psi2, nt, dt, k2, nx, forward, backward);
        else if (order == 4 && type == 'S')
            T4_S(psi, nt, dt, k2, nx, forward, backward);
        else if (order == 6 && type == 'M')
            T6_M(psi, psi1, psi2, psi3, nt, dt, k2, nx, forward, backward);
        else if (order == 6 && type == 'S')
            T6_S(psi, nt, dt, k2, nx, forward, backward);
        else if (order == 8 && type == 'M')
            T8_M(psi, psi1, psi2, psi3, psi4, nt, dt, k2, nx, forward, backward);
        else if (order == 8 && type == 'S')
            T8_S(psi, nt, dt, k2, nx, forward, backward);
        else
        {
            printf("Unknown order. Error.\n"); //Exit somehow
            return 1;
        }
    }

    // Print output
    printf("Printing output.\n");
    char im_f[20];
    char r_f[20];
    char x_f[20];
    sprintf(im_f, "output/psi_i.bin");
    sprintf(r_f, "output/psi_r.bin");
    sprintf(x_f, "output/x.bin");
    
    fp=fopen(im_f, "wb");
    fwrite(psi_i, sizeof(double), nx*nt, fp);
    fclose(fp);
    fp=fopen(r_f, "wb");
    fwrite(psi_r, sizeof(double), nx*nt, fp);
    fclose(fp);
    fp=fopen(x_f, "wb");
    fwrite(x, sizeof(double), nx, fp);
    fclose(fp);
    
    // Clean up
    printf("Cleaninig up.\n");
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(psi);
    free(psi_i); free(psi_r);
    free(x);
    free(k2);
    free(k);

    if (order > 2)
    {
        fftw_free(psi1);
        fftw_free(psi2);
    }
    if (order > 4)
        fftw_free(psi3);
    if (order > 6)
        fftw_free(psi4);

    return 0;
}

