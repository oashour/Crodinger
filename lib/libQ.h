#ifndef LIB_H
#define LIB_H

// Libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <quadmath.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <assert.h>

// Index Macro
#define ind(i,j)  ((i)*sim->nx+(j))		//[i  ,j] 

#define PEREGRINE 1
#define BACKGROUND 2
// PI
#define M_PI 3.14159265358979323846264338327

// Special simulation structure
struct simulation
{
    __float128 dt;
    __float128 tm;
    __float128 l;
    __float128 A1;
    __float128 q;
    int order;
    char type;
    int nx;
    bool print_psi;
    bool print_spectrum;
    int spectrum_sampling;
    int psi_sampling;
    double nu;
    char f_psi_i[256];
    char f_psi_r[256];
    char f_spectrum[256];
    char f_param[256];
    char f_x[256];
    int r_size;
    int s_size;
    __float128 dx;
    int nt;
    __float128 Omega;
    __float128 A0;
    int initial;
};

// IO operations
struct simulation *read_param(void);
void simulation_destroy(struct simulation *sim);
void simulation_print(struct simulation *sim);
void print_param(struct simulation *sim);

// Array operations
void arr_mult(fftwq_complex *psi, __float128 mult, int size);
void arr_add(fftwq_complex *psi, fftwq_complex *psi1, fftwq_complex *psi2, int size);

// Second order
void T2(fftwq_complex *psi, fftwq_complex *psi_o, int nt, __float128 dt, __float128 *k2, int nx,
        fftwq_plan forward, fftwq_plan backward);

// Multi-product integrators
void T4_M(fftwq_complex *psi, fftwq_complex *psi1, fftwq_complex *psi2, int nt, __float128 dt, 
          __float128 *k2, int nx, fftwq_plan forward, fftwq_plan backward);
void T6_M(fftwq_complex *psi, fftwq_complex *psi1, fftwq_complex *psi2, fftwq_complex *psi3,
          int nt, __float128 dt, __float128 *k2, int nx, fftwq_plan forward, fftwq_plan backward);
void T8_M(fftwq_complex *psi, fftwq_complex *psi1, fftwq_complex *psi2, fftwq_complex *psi3,
          fftwq_complex *psi4, int nt, __float128 dt, __float128 *k2, int nx, fftwq_plan forward, 
          fftwq_plan backward);

// Symplectic integrators
void T4_S(fftwq_complex *psi, int nt, __float128 dt, __float128 *k2, int nx, 
          fftwq_plan forward, fftwq_plan backward);
void T6_S(fftwq_complex *psi, int nt, __float128 dt, __float128 *k2, int nx, 
          fftwq_plan forward, fftwq_plan backward);
void T8_S(fftwq_complex *psi, int nt, __float128 dt, __float128 *k2, int nx,
          fftwq_plan forward, fftwq_plan backward);

#endif //LIB_H
