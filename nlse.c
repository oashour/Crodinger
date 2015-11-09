#include "lib.h"
#include <assert.h>// in .h file, why do I need it here? fix.

int main(int argc, char *argv[])
{
    // Read parameters
    struct simulation *sim = read_param();
    assert(sim != NULL);

    printf("Size of MYTYPE=%zd\n", sizeof(MYTYPE));

    printf("\n-----------------------------------\n");
    printf("Running Program.\n");
    printf(  "-----------------------------------\n");
    // Print basic info about simulation to stdout
    simulation_print(sim);

    // Allocate the arrays
    printf("Allocating arrays.\n");
    MYTYPE *psi_r, *psi_i, *spectrum;
    MYFFTW(_complex) *psi_f;
    MYTYPE *x  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx);
    MYTYPE *k2 = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx);
    MYTYPE *k  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx);
    MYFFTW(_complex) *psi = (MYFFTW(_complex)*) MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);
    
    // Full spectrum
    spectrum = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx*sim->print_size);
   
    // Step FFT
    psi_f = (MYFFTW(_complex)*)MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);    
    psi_r = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx*sim->print_size);    // Real part
    psi_i = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->nx*sim->print_size);   // Imaginary part
    
    // Allocating arrays for energy calculation
    MYTYPE *ke  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->print_size);
    MYTYPE *pe  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->print_size);
    MYTYPE *E  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->print_size);
    MYTYPE *dE  = (MYTYPE*)malloc(sizeof(MYTYPE) * sim->print_size);


    // Allocating extra arrays for high orders
    printf("Allocating extra arrays.\n");
    MYFFTW(_complex) *psi1, *psi2, *psi3, *psi4;
    if (sim->order >= 4)                          // Special arrays for order 4, 6, 8 
    {
        psi1 = (MYFFTW(_complex)*)MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);
        psi2 = (MYFFTW(_complex)*)MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);
    }
    if (sim->order >= 6)                              // Special arrays for order 6, 8 
        psi3 = (MYFFTW(_complex)*)MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);
    if (sim->order >= 8)                              // Special arrays for order 8  
        psi4 = (MYFFTW(_complex)*)MYFFTW(_malloc)(sizeof(MYFFTW(_complex)) * sim->nx);

    // Create transform plans
    printf("Creating plans.\n");
    MYFFTW(_plan) forward, backward, f0_plan;
    forward = MYFFTW(_plan_dft_1d)(sim->nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = MYFFTW(_plan_dft_1d)(sim->nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    f0_plan = MYFFTW(_plan_dft_1d)(sim->nx, psi, psi_f, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create wave number. '
    //Flipped to force-shift 0 frequency component to the middle when in k-domain.
    printf("Preparing wavenumber.\n");
    MYTYPE dk=2*M_PI/sim->nx/sim->dx;
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
    MYFFTW(_destroy_plan)(forward);
    MYFFTW(_destroy_plan)(backward);
    MYFFTW(_destroy_plan)(f0_plan);

    // Free base arrays
    MYFFTW(_free)(psi);
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
    MYFFTW(_free)(psi_f);

    // Free high order arrays
    if (sim->order >= 4)
    {
        MYFFTW(_free)(psi1);
        MYFFTW(_free)(psi2);
    }
    if (sim->order >= 6)
        MYFFTW(_free)(psi3);
    if (sim->order >= 8)
        MYFFTW(_free)(psi4);

    // Free simulation struct
    simulation_destroy(sim);

    return 0;
}

