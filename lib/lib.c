#include "lib.h"

void arr_mult(fftw_complex *psi, double mult, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi[i]*mult;
}

void arr_add(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi1[i]+psi2[i];
}

void T2(fftw_complex *psi, fftw_complex *psi_o, int nt, double dt, double *k2, int nx, fftw_plan forward, fftw_plan backward) 
{
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabs(psi[i]) * cabs(psi[i]) * dt/2)*psi[i];
    // Forward transform
    fftw_execute_dft(forward, psi_o, psi_o);
    // Solve linear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(-I * k2[i]/2 * dt)*psi_o[i];
    // Backward transform
    fftw_execute_dft(backward, psi_o, psi_o);
    // Normalize the transform
    arr_mult(psi_o, 1.0/nx, nx);
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabs(psi_o[i]) * cabs(psi_o[i]) * dt/2)*psi_o[i];
}

void T4_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, int nt, double dt, 
          double *k2, int nx, fftw_plan forward, fftw_plan backward)
{
	    T2(psi, psi1, nt, dt/2, k2, nx, forward, backward);
	    T2(psi1, psi1, nt, dt/2, k2, nx, forward, backward);
        arr_mult(psi1, 4.0/3.0, nx);

	    T2(psi, psi2, nt, dt, k2, nx, forward, backward);
        arr_mult(psi2, -1.0/3.0, nx);

        arr_add(psi, psi1, psi2, nx);
}

void T6_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, fftw_complex *psi3,
          int nt, double dt, double *k2, int nx, fftw_plan forward, fftw_plan backward)
{
        T2(psi, psi1, nt, dt/3, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/3, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/3, k2, nx, forward, backward);
        arr_mult(psi1, 81.0/40.0, nx);

        T2(psi, psi2, nt, dt/2, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/2, k2, nx, forward, backward);
        arr_mult(psi2, -16.0/15.0, nx);

        T2(psi, psi3, nt, dt, k2, nx, forward, backward);
        arr_mult(psi3, 1.0/24.0, nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
}

void T8_M(fftw_complex *psi, fftw_complex *psi1, fftw_complex *psi2, fftw_complex *psi3,
          fftw_complex *psi4, int nt, double dt, double *k2, int nx, fftw_plan forward, 
          fftw_plan backward)
{
        T2(psi, psi1, nt, dt/4, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4, k2, nx, forward, backward);
        arr_mult(psi1, 1024.0/315.0, nx);

        T2(psi, psi2, nt, dt/3, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/3, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/3, k2, nx, forward, backward);
        arr_mult(psi2, -729.0/280.0, nx);

        T2(psi, psi3, nt, dt/2, k2, nx, forward, backward);
        T2(psi3, psi3, nt, dt/2, k2, nx, forward, backward);
        arr_mult(psi3, 16.0/45.0, nx);

        T2(psi, psi4, nt, dt, k2, nx, forward, backward);
        arr_mult(psi4, -1.0/360.0, nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
        arr_add(psi, psi, psi4, nx);
}

void T4_S(fftw_complex *psi, int nt, double dt, double *k2, int nx, 
          fftw_plan forward, fftw_plan backward) 
{
        double s = pow(2, 1.0/3.0);
        double os = 1/(2-s);

        double ft = os;
        double bt = -s*os;

        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, bt*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
}

void T6_S(fftw_complex *psi, int nt, double dt, double *k2, int nx, 
          fftw_plan forward, fftw_plan backward) 
{
        double s = pow(2, 1.0/5.0);
        double os = 1/(2-s);

        double ft = os;
        double bt = -s*os;

        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T4_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

void T8_S(fftw_complex *psi, int nt, double dt, double *k2, int nx, 
          fftw_plan forward, fftw_plan backward) 
{
        double s = pow(2, 1.0/7.0);
        double os = 1/(2-s);

        double ft = os;
        double bt = -s*os;

        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T6_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

