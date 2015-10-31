#include "lib.h"

#define PEREGRINE 1
#define BACKGROUND 2

int main(int argc, char *argv[])
{
    // Simulation parameters
    double dt, tm, l; 
    double A1 = 0, q = 0;
    int order, nx, spectrum_sampling, psi_sampling, l_mult, initial; 
    char temp1, temp2, an_choice;
    char type;
    bool print_psi = 0;
    bool print_spectrum = 0;
    double nu = 0;
    char f_psi_i[256], f_psi_r[256], f_spectrum[256], f_param[256], f_x[256];

    // Prepare simulation
    // Basic grid parameters
    printf("dt: ");                         // Grid temporal spacing          
    scanf("%lf", &dt);
    printf("nx: ");                         // Number of Fourier modes
    scanf("%d", &nx);                       
    printf("tm: ");                         // Maximum time
    scanf("%lf", &tm);

    // Initial Conditions
    printf("For the initial wavefunction, [1] denotes "
           "Peregrine soliton, [2] denotes background.\n"
           "Initial wavefunction code: ");
    scanf("%d", &initial);
    if (initial == BACKGROUND)
    {
        printf("A1: ");                         // Cosine amplitude 
        scanf("%lf", &A1);
        printf("Do you want to pick a or n? [a/n]");                          // Breather parameter
        scanf(" %c", &an_choice);
        if (an_choice == 'a')
        {
            printf("a: ");                          // Breather parameter
            scanf("%lf", &q);
        }
        else if (an_choice == 'n')
        {
            printf("n: ");                          // Breather parameter
            scanf("%lf", &nu);
            q = 0.5*(1-1/(nu*nu));
        }
        else
        {
            printf("Unknown option. Exiting\n");
            return 1;
        }

        printf("Length Multiple: ");            // Breather parameter
        scanf("%d", &l_mult);
        l = l_mult*M_PI/sqrt(1-2*q);		    // Spatial period/twice box length
    }
    else if (initial == PEREGRINE)
    {
        printf("L: ");
        scanf("%lf", &l);
    }
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
    print_psi = (temp1 == 'y');
    if (print_psi)
    {
        printf("Result sampling: ");       // If yes, how often
        scanf("%d", &psi_sampling);
        printf("Real psi file path: ");
        scanf("%s" , f_psi_r);
        printf("Imaginary psi file path: ");
        scanf("%s" , f_psi_i);
    }
    
    // Spectrum output 
    printf("Print spectrum? ");            // Print spectrum or not
    scanf(" %c", &temp2);
    print_spectrum = (temp2 == 'y');
    if (print_spectrum)
    {
        printf("Spectrum sampling: ");     // If yes, how often
        scanf("%d", &spectrum_sampling);
        printf("Spectrum file path: ");
        scanf("%s" , f_spectrum);
    }
    if (print_spectrum | print_psi)
    {
        printf("Param file path: ");
        scanf("%s" , f_param);
        printf("x file path: ");
        scanf("%s" , f_x);
    }

    printf("\n-----------------------------------\n");
    printf("Running Program.\n");
    printf(  "-----------------------------------\n");
    // Derived parameters
    int r_size = 0, s_size = 0;
    const int nt = tm/dt;   		   // Number of temporal nodes
    const double dx = (l / nx);		   // Spatial step size
    double Omega, A0;
    if (initial == BACKGROUND)
    {
        Omega = 2*sqrt(1-2*q);         // Fundamental frequency
        A0 = sqrt(1-2*A1*A1);          // Normalization factor
    }
    if (print_psi)
        r_size = nt/psi_sampling;      // Total size of results array
    if (print_spectrum)
        s_size = nt/spectrum_sampling; // Total size of spectrum array
    
    // Print parameter file
    // Print basic info about simulation
    printf("dt:    %.5f\n"
           "nx:    %d\n"
           "tm:    %.3f\n"
           "A1:    %.1e\n"
           "a:     %.8f\n"
           "Order: %d\n"
           "Type:  %c\n"
           "Sn:    %d\n"
           "Pn:    %d\n", dt, nx, tm, A1, q, order, type, s_size, r_size);

    printf("\nPrinting param.txt file.\n");
    FILE *fp = fopen(f_param, "w");
    fprintf(fp, "%.13f\n"
                "%d\n"
                "%.13f\n"
                "%.13f\n"
                "%.13f\n"
                "%d\n"
                "%d\n"
                "%d\n"
                "%c\n", dt, nx, tm, A1, q, s_size, r_size, order, type);
    fclose(fp);

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
    if (order >= 4)                              // Special arrays for order 4, 6, 8 
    {
        psi1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
        psi2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
    }
    if (order >= 6)                              // Special arrays for order 6, 8 
        psi3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx);
    if (order >= 8)                              // Special arrays for order 8  
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
        if (initial == PEREGRINE)
            psi[i] = 1 - (4/(1+4*x[i]*x[i])) + 0*I;
        else if (initial == BACKGROUND)
	        psi[i] = (A0+2*A1*cos(Omega*x[i])) + 0*I;   // Initial WF 
        else
        {
            printf("Unidentified initial WF code. Exiting");
            return 1;
        }

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
                spectrum[ind(h,j)] = cabs(psi_f[j])/nx;   // Save spectrum
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
        fp=fopen(f_psi_i, "wb");
        fwrite(psi_i, sizeof(double), nx*r_size, fp);
        fclose(fp);
        fp=fopen(f_psi_r, "wb");
        fwrite(psi_r, sizeof(double), nx*r_size, fp);
        fclose(fp);
        fp=fopen(f_x, "wb");
        fwrite(x, sizeof(double), nx, fp);
        fclose(fp);
    } 

    if (print_spectrum)
    {
        printf("Printing spectrum.\n");
        fp=fopen(f_spectrum, "wb");
        fwrite(spectrum, sizeof(double), nx*s_size, fp);
        fclose(fp);
    }

    // Clean up
    printf("Cleaning up.\n");
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
    if (order >= 4)
    {
        fftw_free(psi1);
        fftw_free(psi2);
    }
    if (order >= 6)
        fftw_free(psi3);
    if (order >= 8)
        fftw_free(psi4);

    return 0;
}

