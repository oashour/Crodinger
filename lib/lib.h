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

#define MYTYPE double 
#define MYFFTW(func) fftw ## func 
#define D(num) num ## 0 
#define FS "lf"

// PI
#define M_PI 3.14159265358979323846264338327

// Special simulation structure
struct simulation
{
    MYTYPE dt;
    MYTYPE tm;
    MYTYPE l;
    MYTYPE A1;
    MYTYPE q;
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
    MYTYPE dx;
    int nt;
    MYTYPE Omega;
    MYTYPE A0;
    int initial;
};

// IO operations
struct simulation *read_param(void);
void simulation_destroy(struct simulation *sim);
void simulation_print(struct simulation *sim);
void print_output(MYTYPE *psi_i, MYTYPE *psi_r, MYTYPE *x, 
                  MYTYPE *spectrum, MYTYPE *ke, MYTYPE *pe, 
                  MYTYPE *E, MYTYPE *dE, struct simulation *sim);

// Array operations
void arr_mult(MYFFTW(_complex) *psi, MYTYPE mult, int size);
void arr_mult_l(MYTYPE *out, MYTYPE *in1, MYTYPE *in2, int size);
void arr_add(MYFFTW(_complex) *out, MYFFTW(_complex) *in1, MYFFTW(_complex) *in2, int size);
MYTYPE arr_sum_l(MYTYPE *in, int size);

// Second order
void T2(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi_o, int nt, MYTYPE dt, MYTYPE *k2, int nx,
        MYFFTW(_plan) forward, MYFFTW(_plan) backward);

// Multi-product integrators
void T4_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, int nt, MYTYPE dt, 
          MYTYPE *k2, int nx, MYFFTW(_plan) forward, MYFFTW(_plan) backward);
void T6_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, MYFFTW(_complex) *psi3,
          int nt, MYTYPE dt, MYTYPE *k2, int nx, MYFFTW(_plan) forward, MYFFTW(_plan) backward);
void T8_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, MYFFTW(_complex) *psi3,
          MYFFTW(_complex) *psi4, int nt, MYTYPE dt, MYTYPE *k2, int nx, MYFFTW(_plan) forward, 
          MYFFTW(_plan) backward);

// Symplectic integrators
void T4_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx, 
          MYFFTW(_plan) forward, MYFFTW(_plan) backward);
void T6_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx, 
          MYFFTW(_plan) forward, MYFFTW(_plan) backward);
void T8_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx,
          MYFFTW(_plan) forward, MYFFTW(_plan) backward);

// Energy calculation
void energy(MYFFTW(_complex) *psi, MYTYPE *psi_r, MYTYPE *psi_i, MYTYPE *spectrum,
            MYTYPE *ke, MYTYPE *pe, MYTYPE *E, MYTYPE *dE,
            MYFFTW(_plan) f0_plan, MYFFTW(_complex) *psi_f, MYTYPE *k2, int m, int nx);

#endif //LIB_H
