#include "libL.h"

void simulation_destroy(struct simulation *sim)
{
    assert(sim != NULL);

    free(sim);
}

void simulation_print(struct simulation *sim)
{
    printf("Grid parameters: \n");
    printf("Temporal spacing     (dt): %.4Lf\n", sim->dt);
    printf("Spatial spacing      (dx): %.4Lf\n", sim->dx);
    printf("# x nodes            (nx): %d  \n",  sim->nx);
    printf("# t nodes            (nt): %d  \n",  sim->nt);
    printf("Max time             (tm): %.2Lf\n", sim->tm);
    printf("Box size             (l) : %.4Lf\n", sim->l);
    printf("\n");

    printf("Algorithmic information: \n");
    printf("Algorithm:           (--): %d%c\n", sim->order, sim->type);
    printf("\n");
        
    printf("Initial wave function: \n");
    printf("Modulation Amplitude (A1): %.4Lf\n", sim->A1);
    printf("Constant Background  (A0): %.4Lf\n", sim->A0);
    printf("Fundamental Freq     (Om): %.4Lf\n", sim->Omega);
    printf("AB Parameter         (a) : %.4Lf\n", sim->q);
    printf("# unstable modes     (nu): %d  \n", sim->nu? (int)sim->nu:-1);
    printf("\n");
    
    printf("Output information: \n");
    printf("Print spectrum:      (--): %s  \n", sim->print_spectrum? "yes": "no");
    printf("Print psi:           (--): %s  \n", sim->print_psi? "yes": "no");
    printf("Spectrum samp. rate  (ss): %d  \n", sim->spectrum_sampling);
    printf("Psi samp. rate       (ps): %d  \n", sim->psi_sampling);
    printf("Spectrum # t-points  (sn): %d  \n", sim->s_size);
    printf("Psi # t-points       (pn): %d  \n", sim->r_size);
    printf("\n");

    printf("Output files: \n");
    printf("Real psi file:      (--): %s  \n", sim->f_psi_r);
    printf("Imag psi file:      (--): %s  \n", sim->f_psi_i);
    printf("Spectrum file:      (--): %s  \n", sim->f_spectrum);
    printf("Parameters file:    (--): %s  \n", sim->f_param);
    printf("x-array file:       (--): %s  \n", sim->f_x);
}


void print_param(struct simulation *sim)
{
    FILE *fp = fopen(sim->f_param, "w");

    fprintf(fp, "%.13Lf\n", sim->dt);
    fprintf(fp, "%d\n",    sim->nx);
    fprintf(fp, "%.13Lf\n", sim->tm);
    fprintf(fp, "%.13Lf\n", sim->A1);
    fprintf(fp, "%.13Lf\n", sim->q);
    fprintf(fp, "%d\n",    sim->s_size);
    fprintf(fp, "%d\n",    sim->r_size);
    fprintf(fp, "%d\n",    sim->order);
    fprintf(fp, "%c\n",    sim->type);

    fclose(fp);
}

struct simulation *read_param(void)
{
    // Create simulation structure
    struct simulation *param = malloc(sizeof(struct simulation));

    // Temporary valuables to be read into
    long double dt, tm, l; 
    long double A1 = 0, q = 0;
    int order, nx, spectrum_sampling = 0, psi_sampling = 0, l_mult, initial; 
    char temp1, temp2, an_choice;
    char type;
    bool print_psi = 0;
    bool print_spectrum = 0;
    long double nu = 0;
    char f_psi_i[256], f_psi_r[256], f_spectrum[256], f_param[256], f_x[256];
    
    // Basic grid parameters
    printf("dt: ");                         // Grid temporal spacing          
    scanf("%Lf", &dt);
    printf("nx: ");                         // Number of Fourier modes
    scanf("%d", &nx);                       
    printf("tm: ");                         // Maximum time
    scanf("%Lf", &tm);

    // Initial Conditions
    printf("For the initial wavefunction, [1] denotes "
           "Peregrine soliton, [2] denotes background.\n"
           "Initial wavefunction code: ");
    scanf("%d", &initial);
    if (initial == BACKGROUND)
    {
        printf("A1: ");                         // Cosine amplitude 
        scanf("%Lf", &A1);
        printf("Do you want to pick a or n? [a/n]");                          // Breather parameter
        scanf(" %c", &an_choice);
        if (an_choice == 'a')
        {
            printf("a: ");                          // Breather parameter
            scanf("%Lf", &q);
        }
        else if (an_choice == 'n')
        {
            printf("n: ");                          // Breather parameter
            scanf("%Lf", &nu);
            q = 0.5*(1-1/(nu*nu));
        }
        else
        {
            printf("Unknown option. Exiting\n");
            return NULL;
        }

        printf("Length Multiple: ");            // Breather parameter
        scanf("%d", &l_mult);
        l = l_mult*M_PI/sqrt(1-2*q);		    // Spatial period/twice box length
    }
    else if (initial == PEREGRINE)
    {
        printf("L: ");
        scanf("%Lf", &l);
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
    
    // Derived parameters
    int r_size = 0, s_size = 0;
    const int nt = tm/dt;   	// Number of temporal nodes
    const long double dx = (l/nx);		   // Spatial step size
    long double Omega = 0, A0 = 0;
    if (initial == BACKGROUND)
    {
        Omega = 2*sqrt(1-2*q);         // Fundamental frequency
        A0 = sqrt(1-2*A1*A1);          // Normalization factor
    }
    if (print_psi)
        r_size = nt/psi_sampling;      // Total size of results array
    if (print_spectrum)
        s_size = nt/spectrum_sampling; // Total size of spectrum array
    

    param->dt = dt;
    param->tm = tm;
    param->l  = l;
    param->A1 = A1;
    param->q  = q;
    param->order = order;
    param->nx = nx;
    param->type = type;
    param->spectrum_sampling = spectrum_sampling;
    param->psi_sampling = psi_sampling;
    param->print_psi = print_psi;
    param->print_spectrum = print_spectrum;
    param->nu = nu;
    strcpy(param->f_psi_i, f_psi_i); 
    strcpy(param->f_psi_r, f_psi_r); 
    strcpy(param->f_spectrum, f_spectrum); 
    strcpy(param->f_param, f_param); 
    strcpy(param->f_x, f_x); 
    param->r_size = r_size;
    param->s_size = s_size;
    param->dx = dx;
    param->nt = nt;
    param->Omega = Omega;
    param->A0 = A0;
    param->initial = initial;

    return param;
}

void arr_mult(fftwl_complex *psi, long double mult, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi[i]*mult;
}

void arr_add(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi1[i]+psi2[i];
}

void T2(fftwl_complex *psi, fftwl_complex *psi_o, int nt, long double dt, long double *k2, int nx, fftwl_plan forward, fftwl_plan backward) 
{
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabsl(psi[i]) * cabsl(psi[i]) * dt/2.0Q)*psi[i];
    // Forward transform
    fftwl_execute_dft(forward, psi_o, psi_o);
    // Solve linear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(-I * k2[i]/2.0Q * dt)*psi_o[i];
    // Backward transform
    fftwl_execute_dft(backward, psi_o, psi_o);
    // Normalize the transform
    arr_mult(psi_o, 1.0Q/nx, nx);
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabsl(psi_o[i]) * cabsl(psi_o[i]) * dt/2.0Q)*psi_o[i];
}

void T4_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, int nt, long double dt, 
          long double *k2, int nx, fftwl_plan forward, fftwl_plan backward)
{
	    T2(psi, psi1, nt, dt/2.0Q, k2, nx, forward, backward);
	    T2(psi1, psi1, nt, dt/2.0Q, k2, nx, forward, backward);
        arr_mult(psi1, 4.0Q/3.0Q, nx);

	    T2(psi, psi2, nt, dt, k2, nx, forward, backward);
        arr_mult(psi2, -1.0Q/3.0Q, nx);

        arr_add(psi, psi1, psi2, nx);
}

void T6_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, fftwl_complex *psi3,
          int nt, long double dt, long double *k2, int nx, fftwl_plan forward, fftwl_plan backward)
{
        T2(psi, psi1, nt, dt/3.0Q, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/3.0Q, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/3.0Q, k2, nx, forward, backward);
        arr_mult(psi1, 81.0Q/40.0Q, nx);

        T2(psi, psi2, nt, dt/2.0Q, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/2.0Q, k2, nx, forward, backward);
        arr_mult(psi2, -16.0Q/15.0Q, nx);

        T2(psi, psi3, nt, dt, k2, nx, forward, backward);
        arr_mult(psi3, 1.0Q/24.0Q, nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
}

void T8_M(fftwl_complex *psi, fftwl_complex *psi1, fftwl_complex *psi2, fftwl_complex *psi3,
          fftwl_complex *psi4, int nt, long double dt, long double *k2, int nx, fftwl_plan forward, 
          fftwl_plan backward)
{
        T2(psi, psi1, nt, dt/4.0Q, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4.0Q, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4.0Q, k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/4.0Q, k2, nx, forward, backward);
        arr_mult(psi1, 1024.0Q/315.0Q, nx);

        T2(psi, psi2, nt, dt/3.0Q, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/3.0Q, k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/3.0Q, k2, nx, forward, backward);
        arr_mult(psi2, -729.0Q/280.0Q, nx);

        T2(psi, psi3, nt, dt/2.0Q, k2, nx, forward, backward);
        T2(psi3, psi3, nt, dt/2.0Q, k2, nx, forward, backward);
        arr_mult(psi3, 16.0Q/45.0Q, nx);

        T2(psi, psi4, nt, dt, k2, nx, forward, backward);
        arr_mult(psi4, -1.0Q/360.0Q, nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
        arr_add(psi, psi, psi4, nx);
}

void T4_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx, 
          fftwl_plan forward, fftwl_plan backward) 
{
        long double s = pow(2, 1.0Q/3.0Q);
        long double os = 1.0Q/(2.0Q-s);

        long double ft = os;
        long double bt = -s*os;

        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, bt*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
}

void T6_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx, 
          fftwl_plan forward, fftwl_plan backward) 
{
        long double s = pow(2, 1.0Q/5.0Q);
        long double os = 1.0Q/(2.0Q-s);

        long double ft = os;
        long double bt = -s*os;

        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T4_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

void T8_S(fftwl_complex *psi, int nt, long double dt, long double *k2, int nx, 
          fftwl_plan forward, fftwl_plan backward) 
{
        long double s = pow(2, 1.0Q/7.0Q);
        long double os = 1.0Q/(2.0Q-s);

        long double ft = os;
        long double bt = -s*os;

        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T6_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

