#include "libL.h"
#include <assert.h>// in .h file, why do I need it here? fix.

int main(int argc, char *argv[])
{

    // Read parameters
    struct simulation *sim = read_param();
    assert(sim != NULL);

    printf("\n-----------------------------------\n");
    printf("Running Program.\n");
    printf(  "-----------------------------------\n");
    // Print basic info about simulation to stdout
    simulation_print(sim);

    // Allocate the arrays
    printf("Allocating arrays.\n");
    long double *psi_r, *psi_i, *spectrum;
    fftwl_complex *psi_f;
    long double *x  = (long double*)malloc(sizeof(long double) * sim->nx);
    long double *k2 = (long double*)malloc(sizeof(long double) * sim->nx);
    long double *k  = (long double*)malloc(sizeof(long double) * sim->nx);
    fftwl_complex *psi = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * sim->nx);
    
    // Full spectrum
    spectrum = (long double*)malloc(sizeof(long double) * sim->nx*sim->print_size);
   
    // Step FFT
    psi_f = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * sim->nx);    
    psi_r = (long double*)malloc(sizeof(long double) * sim->nx*sim->print_size);    // Real part
    psi_i = (long double*)malloc(sizeof(long double) * sim->nx*sim->print_size);   // Imaginary part
    
    // Allocating arrays for energy calculation
    long double *ke  = (long double*)malloc(sizeof(long double) * sim->print_size);
    long double *pe  = (long double*)malloc(sizeof(long double) * sim->print_size);
    long double *E  = (long double*)malloc(sizeof(long double) * sim->print_size);
    long double *dE  = (long double*)malloc(sizeof(long double) * sim->print_size);


    // Allocating extra arrays for high orders
    printf("Allocating extra arrays.\n");
    fftwl_complex *psi1, *psi2, *psi3, *psi4;
    if (sim->order >= 4)                          // Special arrays for order 4, 6, 8 
    {
        psi1 = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * sim->nx);
        psi2 = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * sim->nx);
    }
    if (sim->order >= 6)                              // Special arrays for order 6, 8 
        psi3 = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * sim->nx);
    if (sim->order >= 8)                              // Special arrays for order 8  
        psi4 = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) * sim->nx);

    // Create transform plans
    printf("Creating plans.\n");
    fftwl_plan forward, backward, f0_plan;
    forward = fftwl_plan_dft_1d(sim->nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftwl_plan_dft_1d(sim->nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    f0_plan = fftwl_plan_dft_1d(sim->nx, psi, psi_f, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create wave number. '
    //Flipped to force-shift 0 frequency component to the middle when in k-domain.
    printf("Preparing wavenumber.\n");
    long double dk=2*M_PI/sim->nx/sim->dx;
    for(int i = sim->nx/2; i >= 0; i--) 
        k[sim->nx/2-i]=(sim->nx/2-i)*dk;
    for(int i = sim->nx/2+1; i < sim->nx; i++)
        k[i]=(i-sim->nx)*dk; 

    // Initial conditions
    printf("Setting up initial conditions.\n");
    for (int i = 0; i < sim->nx; i++)
    {
        x[i] = (i-sim->nx/2)*sim->dx;               // x array
        k2[i] = k[i]*k[i];                          // Square wave number
        if (sim->initial == PEREGRINE)
            psi[i] = 1 - (4/(1+4*x[i]*x[i])) + 0*I;
        else if (sim->initial == BACKGROUND)
            psi[i] = (sim->A0+2*sim->A1*cos(sim->Omega*x[i])) + 0*I;   // Initial WF 
        else
        {
            printf("Unidentified initial WF code. Exiting");
            return 1;
        }

    }
   
    // Start time evolution
    int m = 0;                                      // Spectrum/psi size 
    printf("Evolving time.\n");
    for (int i = 0; i < sim->nt; i++)
    {
        // Save results for printing
        if (i % sim->sampling == 0)
        {
            energy(psi, psi_r, psi_i, spectrum, ke, pe, E, dE, f0_plan, psi_f, k2, m, sim->nx);
            m++;
        }

        // Evolve one step in time, depending on algorithm
         if (sim->order == 2)
            T2(psi, psi, sim->nt, sim->dt, k2, sim->nx, forward, backward); 
        else if (sim->order == 4 && sim->type == 'M')
            T4_M(psi, psi1, psi2, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else if (sim->order == 4 && sim->type == 'S')
            T4_S(psi, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else if (sim->order == 6 && sim->type == 'M')
            T6_M(psi, psi1, psi2, psi3, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else if (sim->order == 6 && sim->type == 'S')
            T6_S(psi, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else if (sim->order == 8 && sim->type == 'M')
            T8_M(psi, psi1, psi2, psi3, psi4, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else if (sim->order == 8 && sim->type == 'S')
            T8_S(psi, sim->nt, sim->dt, k2, sim->nx, forward, backward);
        else
        {
            printf("Unknown order. Error.\n");
            return 1;
        }
    }
    

    // Print output
    printf("Printing output.\n");
    print_output(psi_i, psi_r, x, spectrum, ke, pe, E, dE, sim);
    
    // Clean up
    printf("Cleaning up.\n");
    // Destroy plans
    fftwl_destroy_plan(forward);
    fftwl_destroy_plan(backward);
    fftwl_destroy_plan(f0_plan);

    // Free base arrays
    fftwl_free(psi);
    free(x);
    free(k2);
    free(k);
    
    // Free output arrays
    free(spectrum);
    free(pe);
    free(ke);
    free(E);
    free(dE);
    free(psi_i); 
    free(psi_r);
    fftwl_free(psi_f);

    // Free high order arrays
    if (sim->order >= 4)
    {
        fftwl_free(psi1);
        fftwl_free(psi2);
    }
    if (sim->order >= 6)
        fftwl_free(psi3);
    if (sim->order >= 8)
        fftwl_free(psi4);

    // Free simulation struct
    simulation_destroy(sim);

    return 0;
}

