/*s********************************************************************************
* Numerical Solution for the Cubic Nonlinear Schrodinger Equation in (1+1)D	 	  *
* using symmetric split step Fourier method		                              	  *
* Coded by: Omar Ashour, Texas A&M University at Qatar, February 2015.    	      *
* ********************************************************************************/
#include "lib.h"
#define ENERGY 0 

int main(int argc, char *argv[])
{
    const int print_results = 1;
    const int fourier_0 = 0;
    const int n_print = 100;

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

    const int nt = tm/dt;   			        // number of temporal nodes
    const double l = M_PI/sqrt(1-2*q);		    // Spatial Period
    const double dx = (l / nx);			        // spatial step size
    const double Omega = 2*sqrt(1-2*q);
    const double A0 = sqrt(1-2*A1*A1);
    const int size = nt/n_print;
    
    char paramf[20];
    sprintf(paramf, "output/param.txt");
    FILE *fp = fopen(paramf, "w");
    fprintf(fp, "%.13f\n"
                "%d\n"
                "%.13f\n"
                "%.13f\n"
                "%.13f\n"
                "%d\n"
                "%c\n"
                "%d\n", dt, nx, tm, A1, q, order, type, size);
    fclose(fp);

    // Print basic info about simulation
    printf("dt=%f\nnx=%d\ntm=%f\norder=%d\ntype=%c\nq=%f\n", 
            dt, nx, tm, order, type, q);

    // Allocate the arrays
    double *psi_r, *psi_i, *psi_f_0, *psi_f_1, *psi_f_2, *psi_f_3, *psi_f_4, *psi_f_5;
    fftw_complex *psi_f;
    double *x = (double*)malloc(sizeof(double) * nx);
    double *k2 = (double*)malloc(sizeof(double) * nx);
    fftw_complex *psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    if (fourier_0)
    {
        psi_f_0 = (double*) malloc(sizeof(double) * size);
        psi_f_1 = (double*) malloc(sizeof(double) * size);
        psi_f_2 = (double*) malloc(sizeof(double) * size);
        psi_f_3 = (double*) malloc(sizeof(double) * size);
        psi_f_4 = (double*) malloc(sizeof(double) * size);
        psi_f_5 = (double*) malloc(sizeof(double) * size);
        psi_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    }
    if (print_results)
    {
        psi_r = (double*) malloc(sizeof(double) * nx * nt);
        psi_i = (double*) malloc(sizeof(double) * nx * nt);
    }

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
    fftw_plan forward, backward, f0_plan;
    forward = fftw_plan_dft_1d(nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_1d(nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (fourier_0)
        f0_plan = fftw_plan_dft_1d(nx, psi, psi_f, FFTW_FORWARD, FFTW_ESTIMATE);

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
	    psi[i] = (A0+2*A1*cos(Omega*x[i])) + 0*I; 
    }
   
    // Start time evolution
    int h = 0;
    printf("Evolving time.\n");
    for (int i = 0; i < nt; i++)
    {
        if (fourier_0 && (i % n_print == 0))
        {
            fftw_execute(f0_plan);
            psi_f_0[h] = cabs(psi_f[0])/nx;
            psi_f_1[h] = cabs(psi_f[1])/nx;
            psi_f_2[h] = cabs(psi_f[2])/nx;
            psi_f_3[h] = cabs(psi_f[3])/nx;
            psi_f_4[h] = cabs(psi_f[4])/nx;
            psi_f_5[h] = cabs(psi_f[5])/nx;
            h++;
        }
     
        if (print_results)
        {
            for (int j = 0; j < nx; j++)
            {
                psi_i[ind(i,j)] = cimag(psi[j]);
                psi_r[ind(i,j)] = creal(psi[j]);
            }	
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
    printf("%d/%d", size, h);
    // Print output
    if (print_results)
    {
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
    } 

    if (fourier_0)
    {
        fp=fopen("output/psi_f_0.bin", "wb");
        fwrite(psi_f_0, sizeof(double), size, fp);
        fclose(fp);
        fp=fopen("output/psi_f_1.bin", "wb");
        fwrite(psi_f_1, sizeof(double), size, fp);
        fclose(fp);
        fp=fopen("output/psi_f_2.bin", "wb");
        fwrite(psi_f_2, sizeof(double), size, fp);
        fclose(fp);
        fp=fopen("output/psi_f_3.bin", "wb");
        fwrite(psi_f_3, sizeof(double), size, fp);
        fclose(fp);
        fp=fopen("output/psi_f_4.bin", "wb");
        fwrite(psi_f_4, sizeof(double), size, fp);
        fclose(fp);
        fp=fopen("output/psi_f_5.bin", "wb");
        fwrite(psi_f_5, sizeof(double), size, fp);
        fclose(fp);
    }

    // Clean up
    printf("Cleaninig up.\n");
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    if (fourier_0)
        fftw_destroy_plan(f0_plan);
    fftw_free(psi);
    if (fourier_0)
    {
        free(psi_f_0);
        free(psi_f_1);
        free(psi_f_2);
        free(psi_f_3);
        free(psi_f_4);
        free(psi_f_5);
        fftw_free(psi_f);
    }
    if (print_results)
    {
        free(psi_i); 
        free(psi_r);
    }
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

