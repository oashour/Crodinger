#include "libQ.h"
#include <assert.h>// in .h file, why do I need it here? fix.

int main(int argc, char *argv[])
{

    // Read parameters
    struct simulation *sim = read_param();
    assert(sim != NULL);

    printf("\n-----------------------------------\n");
    printf("Running Program.\n");
    printf(  "-----------------------------------\n");
    // Print basic info about simulation
    simulation_print(sim);
    // Print parameter file
     printf("\nPrinting param.txt file.\n");
    print_param(sim);

    // Allocate the arrays
    printf("Allocating arrays.\n");
    __float128 *psi_r, *psi_i, *spectrum;
    fftwq_complex *psi_f;
    __float128 *x  = (__float128*)malloc(sizeof(__float128) * sim->nx);
    __float128 *k2 = (__float128*)malloc(sizeof(__float128) * sim->nx);
    __float128 *k  = (__float128*)malloc(sizeof(__float128) * sim->nx);
    fftwq_complex *psi = (fftwq_complex*) fftwq_malloc(sizeof(fftwq_complex) * sim->nx);
    if (sim->print_spectrum)
    { 
        // Full spectrum
        spectrum = (__float128*)malloc(sizeof(__float128) * sim->nx*sim->s_size);
        // Step FFT
        psi_f = (fftwq_complex*)fftwq_malloc(sizeof(fftwq_complex) * sim->nx);    
    }
    if (sim->print_psi)
    {
        psi_r = (__float128*)malloc(sizeof(__float128) * sim->nx*sim->r_size);    // Real part
        psi_i = (__float128*)malloc(sizeof(__float128) * sim->nx*sim->r_size);   // Imaginary part
    }

    // Allocating extra arrays for high orders
    printf("Allocating extra arrays.\n");
    fftwq_complex *psi1, *psi2, *psi3, *psi4;
    if (sim->order >= 4)                          // Special arrays for order 4, 6, 8 
    {
        psi1 = (fftwq_complex*)fftwq_malloc(sizeof(fftwq_complex) * sim->nx);
        psi2 = (fftwq_complex*)fftwq_malloc(sizeof(fftwq_complex) * sim->nx);
    }
    if (sim->order >= 6)                              // Special arrays for order 6, 8 
        psi3 = (fftwq_complex*)fftwq_malloc(sizeof(fftwq_complex) * sim->nx);
    if (sim->order >= 8)                              // Special arrays for order 8  
        psi4 = (fftwq_complex*)fftwq_malloc(sizeof(fftwq_complex) * sim->nx);

    // Create transform plans
    printf("Creating plans.\n");
    fftwq_plan forward, backward, f0_plan;
    forward = fftwq_plan_dft_1d(sim->nx, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftwq_plan_dft_1d(sim->nx, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (sim->print_spectrum)
        f0_plan = fftwq_plan_dft_1d(sim->nx, psi, psi_f, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create wave number. '
    //Flipped to force-shift 0 frequency component to the middle when in k-domain.
    printf("Preparing wavenumber.\n");
    __float128 dk=2*M_PIq/sim->nx/sim->dx;
    for(int i = sim->nx/2; i >= 0; i--) 
        k[sim->nx/2-i]=(sim->nx/2-i)*dk;
    for(int i = sim->nx/2+1; i < sim->nx; i++)
        k[i]=(i-sim->nx)*dk; 

    // Initial conditions
    printf("Setting up initial conditions.\n");
    for (int i = 0; i < sim->nx; i++)
    {
        x[i] = (i-sim->nx/2)*sim->dx;                         // x array
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
    int h = 0, m = 0;                                      // Spectrum/psi size 
    printf("Evolving time.\n");
    for (int i = 0; i < sim->nt; i++)
    {
        if (sim->print_spectrum && (i % sim->spectrum_sampling == 0))
        {
            fftwq_execute(f0_plan);                         // Get step's spectrum 
            for (int j = 0; j < sim->nx; j++)
                spectrum[ind(h,j)] = cabsq(psi_f[j])/sim->nx;   // Save spectrum
            h++;
        }
     
        // Save results for printing
        if (sim->print_psi && (i % sim->psi_sampling == 0))
        {
            for (int j = 0; j < sim->nx; j++)
            {
                psi_r[ind(m,j)] = crealq(psi[j]);
                psi_i[ind(m,j)] = cimagq(psi[j]);
            }	
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
    double *psi_rd, *psi_id, *x_d;
    FILE *fp;
    if (sim->print_psi)
    {
        psi_rd = (double*)malloc(sizeof(double) * sim->nx*sim->r_size);  
        psi_id = (double*)malloc(sizeof(double) * sim->nx*sim->r_size);   
        x_d    = (double*)malloc(sizeof(double) * sim->nx);    
        
        for (int i = 0; i < sim->nx*sim->r_size; i++)
        {
            psi_id[i] = (double) psi_i[i];
            psi_rd[i] = (double) psi_r[i];
        }

        for (int i = 0; i < sim->nx; i++)
            x_d[i] = (double) x[i];
       
        printf("Printing output.\n");
        fp=fopen(sim->f_psi_i, "wb");
        fwrite(psi_id, sizeof(double), sim->nx*sim->r_size, fp);
        fclose(fp);
        fp=fopen(sim->f_psi_r, "wb");
        fwrite(psi_rd, sizeof(double), sim->nx*sim->r_size, fp);
        fclose(fp);
        fp=fopen(sim->f_x, "wb");
        fwrite(x_d, sizeof(double), sim->nx, fp);
        fclose(fp);
    } 

    // if (print_spectrum)
    // {
        // printf("Printing spectrum.\n");
        // fp=fopen(f_spectrum, "wb");
        // fwrite(spectrum, sizeof(__float128), nx*s_size, fp);
        // fclose(fp);
    // }

    // Clean up
    printf("Cleaning up.\n");
    // Destroy plans
    fftwq_destroy_plan(forward);
    fftwq_destroy_plan(backward);
    if (sim->print_spectrum)
        fftwq_destroy_plan(f0_plan);

    // Free base arrays
    fftwq_free(psi);
    free(x);
    free(k2);
    free(k);
    
    // Free conditional arrays
    if (sim->print_spectrum)
    {
        free(spectrum);
        fftwq_free(psi_f);
    }
    if (sim->print_psi)
    {
        free(psi_i); 
        free(psi_r);
    }
    if (sim->order >= 4)
    {
        fftwq_free(psi1);
        fftwq_free(psi2);
    }
    if (sim->order >= 6)
        fftwq_free(psi3);
    if (sim->order >= 8)
        fftwq_free(psi4);

    simulation_destroy(sim);

    return 0;
}

