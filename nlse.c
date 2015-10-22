#include "lib.h"

int main(int argc, char *argv[])
{
    // Simulation parameters
    double dt, A1, q, tm; 
    int order, nx, spectrum_sampling, psi_sampling; 
    char temp1, temp2;
    char type;

    // Prepare simulation
    // Basic grid parameters
    printf("dt: ");                         // Grid temporal spacing          
    scanf("%lf", &dt);
    printf("tm: ");                         // Maximum time
    scanf("%lf", &tm);
    printf("nx: ");                         // Number of Fourier modes
    scanf("%d", &nx);                       
    printf("A1: ");                         // Cosine amplitude 
    scanf("%lf", &A1);
    printf("a: ");                          // Breather parameter
    scanf("%lf", &q);

    // Determine algorithm 
    printf("Algorithm order: ");            // Pick algorithm
    scanf("%d", &order);
    if(order != 2)
    {
        printf("Integrator type: ");        // Type of integrator
        scanf(" %c", &type);
    }
    else
        type = '-';

    // Psi output 
    printf("Print psi? ");                  // Print psi or not
    scanf(" %c", &temp1);
    const bool print_psi = (temp1 != 'n');
    if (print_psi)
    {
        printf("Result sampling: ");       // If yes, how often
        scanf("%d", &psi_sampling);
    }
    
    // Spectrum output 
    printf("Print spectrum? ");            // Print spectrum or not
    scanf(" %c", &temp2);
    const bool print_spectrum = (temp2 != 'n');
    if (print_spectrum)
    {
        printf("Spectrum sampling: ");     // If yes, how often
        scanf("%d", &spectrum_sampling);
    }

    printf("-----------------------------------.\n");
    printf("Running Program.\n");
    printf("-----------------------------------.\n");
    // Derived parameters
    int r_size, s_size;
    const int nt = tm/dt;   			        // Number of temporal nodes
    const double l = M_PI/sqrt(1-2*q);		    // Spatial period/twice box length
    const double dx = (l / nx);			        // Spatial step size
    const double Omega = 2*sqrt(1-2*q);         // Fundamental frequency
    const double A0 = sqrt(1-2*A1*A1);          // Normalization factor
    if (print_psi)
        r_size = nt/psi_sampling;           // Total size of results array
    if (print_spectrum)
        s_size = nt/spectrum_sampling;           // Total size of spectrum array
    
    // Print parameter file
    printf("Printing param.txt file.\n");
    FILE *fp = fopen("output/param.txt", "w");
    fprintf(fp, "%.13f\n"
                "%d\n"
                "%.13f\n"
                "%.13f\n"
                "%.13f\n"
                "%d\n"
                "%c\n", dt, nx, tm, A1, q, order, type);
    fclose(fp);

    // Print basic info about simulation
    printf("dt:    %.5f\n"
           "nx:    %d\n"
           "tm:    %.3f\n"
           "A1:    %.1e\n"
           "a:     %.8f\n"
           "Order: %d\n"
           "Type:  %c\n", dt, nx, tm, A1, q, order, type);

    // Allocate the arrays
    printf("Allocating arrays.\n");
    double *psi_r, *psi_i, *spectrum;
    fftw_complex *psi_f;
    double *x = (double*)malloc(sizeof(double) * nx);
    double *k2 = (double*)malloc(sizeof(double) * nx);
    double *k = (double*)malloc(sizeof(double) * nx);
    fftw_complex *psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    if (print_spectrum)
    {
        spectrum = (double*)malloc(sizeof(double) * nx*s_size);  // Full spectrum
        psi_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx); // Step FFT
    }
    if (print_psi)
    {
        psi_r = (double*)malloc(sizeof(double) * nx*r_size);     // Real part
        psi_i = (double*)malloc(sizeof(double) * nx*r_size);     // Imaginary part
    }

    // Allocating extra arrays for high orders
    printf("Allocating extra arrays.\n");
    fftw_complex *psi1, *psi2, *psi3, *psi4;
    if (order > 2)                              // Special arrays for order 4, 6, 8 
    {
        psi1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
        psi2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
    }
    if (order > 4)                              // Special arrays for order 6, 8 
        psi3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
    if (order > 6)                              // Special arrays for order 8  
        psi4 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);

    // Create transform plans
    printf("Creating plans.\n");
    fftw_plan forward, backward, f0_plan;
    forward = fftw_plan_dft_1d(nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_1d(nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (print_spectrum)
        f0_plan = fftw_plan_dft_1d(nx, psi, psi_f, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create wave number. '
    //Flipped to force-shift 0 frequency component to the middle when in k-domain.
    printf("Preparing wavenumber.\n");
    double dk=2*M_PI/nx/dx;
    for(int i = nx/2; i >= 0; i--) 
        k[nx/2-i]=(nx/2-i)*dk;
    for(int i = nx/2+1; i < nx; i++)
	    k[i]=(i-nx)*dk; 

    // Initial conditions
    printf("Setting up initial conditions.\n");
    for (int i = 0; i < nx; i++)
    {
        x[i] = (i-nx/2)*dx;                         // x array
	    k2[i] = k[i]*k[i];                          // Square wave number
	    psi[i] = (A0+2*A1*cos(Omega*x[i])) + 0*I;   // Initial WF 
    }
   
    // Start time evolution
    int h = 0, m = 0;                                      // Spectrum/psi size 
    printf("Evolving time.\n");
    for (int i = 0; i < nt; i++)
    {
        if (print_spectrum && (i % spectrum_sampling == 0))
        {
            fftw_execute(f0_plan);                         // Get step's spectrum 
            for (int j = 0; j < nx; j++)
            {
                spectrum[ind(h,j)] = cabs(psi_f[j])/nx;   // Save spectrum
                printf("h/j/ind: %d/%d/%d\n", h, j, ind(h, j));
            }
            h++;
        }
     
        // Save results for printing
        if (print_psi && (i % psi_sampling == 0))
        {
            for (int j = 0; j < nx; j++)
            {
                psi_r[ind(m,j)] = creal(psi[j]);
                psi_i[ind(m,j)] = cimag(psi[j]);
            }	
            m++;
        }

        // Evolve one step in time, depending on algorithm
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
            printf("Unknown order. Error.\n");
            return 1;
        }
    }
    
    // Print output
    if (print_psi)
    {
        printf("Printing output.\n");
        fp=fopen("output/psi_i.bin", "wb");
        fwrite(psi_i, sizeof(double), nx*r_size, fp);
        fclose(fp);
        fp=fopen("output/psi_r.bin", "wb");
        fwrite(psi_r, sizeof(double), nx*r_size, fp);
        fclose(fp);
        fp=fopen("output/x.bin", "wb");
        fwrite(x, sizeof(double), nx, fp);
        fclose(fp);
    } 

    if (print_spectrum)
    {
        printf("Printing spectrum.\n");
        fp=fopen("output/spectrum.bin", "wb");
        fwrite(spectrum, sizeof(double), nt*s_size, fp);
        fclose(fp);
    }

    // Clean up
    printf("Cleaninig up.\n");
    // Destroy plans
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    if (print_spectrum)
        fftw_destroy_plan(f0_plan);

    // Free base arrays
    fftw_free(psi);
    free(x);
    free(k2);
    free(k);
    
    // Free conditional arrays
    if (print_spectrum)
    {
        free(spectrum);
        fftw_free(psi_f);
    }
    if (print_psi)
    {
        free(psi_i); 
        free(psi_r);
    }
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

