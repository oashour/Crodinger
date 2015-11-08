#ifndef LIB_H
#define LIB_H

// Libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <assert.h>

// Index Macro
#define ind(i,j)  ((i)*sim->nx+(j))		//[i  ,j] 
#define ind2(i,j) ((i)*nx+(j))		//[i  ,j] 

#define PEREGRINE 1
#define BACKGROUND 2

// PI
#define M_PI 3.14159265358979323846264338327

// Special simulation structure
struct simulation
{
    long double dt;
    long double tm;
    long double l;
    long double A1;
    long double q;
    int order;
    char type;
    int nx;
    int sampling;
    double nu;
    char f_psi_i[256];
    char f_psi_r[256];
    char f_spectrum[256];
    char f_param[256];
    char f_x[256];
    char f_ke[256];
    char f_pe[256];
    char f_E[256];
    char f_dE[256];
    int print_size;
    long double dx;
    int nt;
    long double Omega;
    long double A0;
    int initial;
};

// IO operations
struct simulation *read_param(void);
void simulation_destroy(struct simulation *sim);
void simulation_print(struct simulation *sim);
void print_output(long double *psi_i, long double *psi_r, long double *x, 
                  long double *spectrum, long double *ke, long double *pe, 
                  long double *E, long double *dE, struct simulation *sim);

// Array operations
void arr_mult(fftwl_complex *psi, long double mult, int size);
void arr_mult_l(long double *out, long double *in1, long double *in2, int size);
void arr_add(fftwl_complex *out, fftwl_complex *in1, fftwl_complex *in2, int size);
long double arr_sum_l(long double *in, int size);

// Second order
void T2(fftwl_complex *psi, fftwl_complex *psi_o, int nt, long double dt, long double *k2, int nx,
        fftwl_plan forward, fftwl_plan backward);

// Multi-product integrators
void T4_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, int nt, long double dt, 
          long double *k2, int nx, fftwl_plan forward, fftwl_plan backward);
void T6_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, fftwl_complex *psi3,
          int nt, long double dt, long double *k2, int nx, fftwl_plan forward, fftwl_plan backward);
void T8_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, fftwl_complex *psi3,
          fftwl_complex *psi4, int nt, long double dt, long double *k2, int nx, fftwl_plan forward, 
          fftwl_plan backward);

// Symplectic integrators
void T4_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx, 
          fftwl_plan forward, fftwl_plan backward);
void T6_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx, 
          fftwl_plan forward, fftwl_plan backward);
void T8_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx,
          fftwl_plan forward, fftwl_plan backward);

// Energy calculation
void energy(fftwl_complex *psi, long double *psi_r, long double *psi_i, long double *spectrum,
            long double *ke, long double *pe, long double *E, long double *dE,
            fftwl_plan f0_plan, fftwl_complex *psi_f, long double *k2, int m, int nx);

#endif //LIB_H
