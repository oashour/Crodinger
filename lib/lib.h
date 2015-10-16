#ifndef LIB_H
#define LIB_H

// Libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

// PI
#define M_PI 3.14159265358979323846264338327

// Index Macro
#define ind(i,j)  ((i)*nx+(j))		//[i  ,j] 

// Array operations
void arr_mult(fftw_complex *psi, double mult, int size);
void arr_add(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, int size);

// Second order
void T2(fftw_complex *psi, fftw_complex *psi_o, int nt, double dt, double *k2, int nx,
        fftw_plan forward, fftw_plan backward);

// Multi-product integrators
void T4_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, int nt, double dt, 
          double *k2, int nx, fftw_plan forward, fftw_plan backward);
void T6_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, fftw_complex *psi3,
          int nt, double dt, double *k2, int nx, fftw_plan forward, fftw_plan backward);
void T8_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, fftw_complex *psi3,
          fftw_complex *psi4, int nt, double dt, double *k2, int nx, fftw_plan forward, 
          fftw_plan backward);

// Symplectic integrators
void T4_S(fftw_complex *psi, int nt, double dt, double *k2, int nx, 
          fftw_plan forward, fftw_plan backward);
void T6_S(fftw_complex *psi, int nt, double dt, double *k2, int nx, 
          fftw_plan forward, fftw_plan backward);
void T8_S(fftw_complex *psi, int nt, double dt, double *k2, int nx,
          fftw_plan forward, fftw_plan backward);

#endif //TIMERS_H
