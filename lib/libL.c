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
    printf("Samp. rate       (ps): %d  \n", sim->sampling);
    printf("# Sampled points (pn): %d  \n", sim->print_size);
    printf("\n");

    printf("Output files: \n");
    printf("Real psi file:       (--): %s  \n", sim->f_psi_r);
    printf("Imag psi file:       (--): %s  \n", sim->f_psi_i);
    printf("Spectrum file:       (--): %s  \n", sim->f_spectrum);
    printf("Parameters file:     (--): %s  \n", sim->f_param);
    printf("KE file:             (--): %s  \n", sim->f_ke);
    printf("PE file              (--): %s  \n", sim->f_pe);
    printf("E file:              (--): %s  \n", sim->f_E);
    printf("dE file:             (--): %s  \n", sim->f_dE);
    printf("x-array file:        (--): %s  \n", sim->f_x);
}

void print_output(long double *psi_i, long double *psi_r, long double *x, 
                  long double *spectrum, long double *ke, long double *pe, 
                  long double *E, long double *dE, struct simulation *sim)
{
    FILE *fp = fopen(sim->f_param, "w");

    fprintf(fp, "%.13Lf\n", sim->dt);
    fprintf(fp, "%d\n",    sim->nx);
    fprintf(fp, "%.13Lf\n", sim->tm);
    fprintf(fp, "%.13Lf\n", sim->A1);
    fprintf(fp, "%.13Lf\n", sim->q);
    fprintf(fp, "%d\n",    sim->print_size);
    fprintf(fp, "%d\n",    sim->order);
    fprintf(fp, "%c\n",    sim->type);

    double *psi_rd, *psi_id, *x_d, *pe_d;
    double *spectrum_d, *ke_d, *dE_d, *E_d;

    psi_rd = (double*)malloc(sizeof(double) * sim->nx*sim->print_size);  
    psi_id = (double*)malloc(sizeof(double) * sim->nx*sim->print_size);   
    x_d    = (double*)malloc(sizeof(double) * sim->nx);    
    pe_d   = (double*)malloc(sizeof(double) * sim->print_size);    
    spectrum_d = (double*)malloc(sizeof(double) * sim->nx*sim->print_size);  
    ke_d = (double*)malloc(sizeof(double) * sim->print_size);  
    E_d = (double*)malloc(sizeof(double) * sim->print_size);  
    dE_d = (double*)malloc(sizeof(double) * sim->print_size);  
    
    for (int i = 0; i < sim->nx*sim->print_size; i++)
    {
        psi_id[i] = (double) psi_i[i];
        psi_rd[i] = (double) psi_r[i];
    }

    for (int i = 0; i < sim->nx; i++)
        x_d[i] = (double) x[i];

    for (int i = 0; i < sim->print_size; i++)
    {
        pe_d[i] = (double) pe[i];
        dE_d[i] = (double) dE[i];
        E_d[i] = (double) E[i];
        ke_d[i] = (double) ke[i];
    }

    for (int i = 0; i < sim->nx*sim->print_size; i++)
        spectrum_d[i] = (double) spectrum[i];

    fp=fopen(sim->f_psi_i, "wb");
    fwrite(psi_id, sizeof(double), sim->nx*sim->print_size, fp);
    fclose(fp);
    fp=fopen(sim->f_psi_r, "wb");
    fwrite(psi_rd, sizeof(double), sim->nx*sim->print_size, fp);
    fclose(fp);
    fp=fopen(sim->f_x, "wb");
    fwrite(x_d, sizeof(double), sim->nx, fp);
    fp=fopen(sim->f_spectrum, "wb");
    fwrite(spectrum_d, sizeof(double), sim->nx*sim->print_size, fp);
    fp=fopen(sim->f_pe, "wb");
    fwrite(pe_d, sizeof(double), sim->print_size, fp);
    fp=fopen(sim->f_ke, "wb");
    fwrite(ke_d, sizeof(double), sim->print_size, fp);
    fp=fopen(sim->f_E, "wb");
    fwrite(E_d, sizeof(double), sim->print_size, fp);
    fp=fopen(sim->f_dE, "wb");
    fwrite(dE_d, sizeof(double), sim->print_size, fp);
    fclose(fp);

    free(ke_d); free(pe_d); free(dE_d); free(E_d);
    free(psi_id); free(spectrum_d); 
    free(psi_rd);
}
struct simulation *read_param(void)
{
    // Create simulation structure
    struct simulation *param = malloc(sizeof(struct simulation));

    // Temporary valuables to be read into
    long double dt, tm, l; 
    long double A1 = 0, q = 0;
    int order, nx, sampling, l_mult, initial; 
    char an_choice;
    char type;
    long double nu = 0;
    char f_psi_i[256], f_psi_r[256], f_spectrum[256], f_param[256], f_x[256];
    char f_pe[256], f_ke[256], f_E[256], f_dE[256];
    
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

    // Output 
    printf("Result sampling: ");       
    scanf("%d", &sampling);
    printf("Real psi file path: ");
    scanf("%s" , f_psi_r);
    printf("Imaginary psi file path: ");
    scanf("%s" , f_psi_i);
    printf("Spectrum file path: ");
    scanf("%s" , f_spectrum);
    printf("Param file path: ");
    scanf("%s" , f_param);
    printf("x file path: ");
    scanf("%s" , f_x);
    printf("PE file path: ");
    scanf("%s" , f_pe);
    printf("KE file path: ");
    scanf("%s" , f_ke);
    printf("E file path: ");
    scanf("%s" , f_E);
    printf("dE file path: ");
    scanf("%s" , f_dE);
    
    // Derived parameters
    int print_size;
    const int nt = tm/dt;   	           // Number of temporal nodes
    const long double dx = (l/nx);		   // Spatial step size
    long double Omega = 0, A0 = 0;
    if (initial == BACKGROUND)
    {
        Omega = 2*sqrt(1-2*q);         // Fundamental frequency
        A0 = sqrt(1-2*A1*A1);          // Normalization factor
    }
    print_size = nt/sampling;      // Total size of results array
    

    param->dt = dt;
    param->tm = tm;
    param->l  = l;
    param->A1 = A1;
    param->q  = q;
    param->order = order;
    param->nx = nx;
    param->type = type;
    param->sampling = sampling;
    param->nu = nu;
    strcpy(param->f_psi_i, f_psi_i); 
    strcpy(param->f_psi_r, f_psi_r); 
    strcpy(param->f_spectrum, f_spectrum); 
    strcpy(param->f_param, f_param); 
    strcpy(param->f_x, f_x); 
    strcpy(param->f_ke, f_ke); 
    strcpy(param->f_pe, f_pe); 
    strcpy(param->f_E, f_E); 
    strcpy(param->f_dE, f_dE); 
    param->print_size = print_size;
    param->dx = dx;
    param->nt = nt;
    param->Omega = Omega;
    param->A0 = A0;
    param->initial = initial;

    return param;
}

void arr_mult_c(fftwl_complex *out, fftwl_complex *in1, fftwl_complex *in2, int size)
{
    for(int i =0; i < size; i++)
        out[i] = in1[i]*in2[i];
}

void arr_mult_l(long double *out, long double *in1, long double *in2, int size)
{
    for(int i =0; i < size; i++)
        out[i] = in1[i]*in2[i];
}

void arr_mult(fftwl_complex *psi, long double mult, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi[i]*mult;
}

long double arr_sum_l(long double *in, int size)
{
    long double out = 0;
    for (int i =0; i < size; i++)
        out += in[i];

    return out;
}

void arr_add_l(long double *out, long double *in1, long double *in2, int size)
{
    for (int i = 0; i < size; i++)
        out[i] = in1[i]+in2[i];
}

void arr_add(fftwl_complex *out, fftwl_complex *in1, fftwl_complex *in2, int size)
{
    for (int i = 0; i < size; i++)
        out[i] = in1[i]+in2[i];
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

void energy(fftwl_complex *psi, long double *psi_r, long double *psi_i, long double *spectrum,
            long double *ke, long double *pe, long double *E, long double *dE,
            fftwl_plan f0_plan, fftwl_complex *psi_f, long double *k2, int m, int nx)
{
    long double *temp1 = (long double*)malloc(sizeof(long double) * nx);
    long double *temp2 = (long double*)malloc(sizeof(long double) * nx);

    for (int j = 0; j < nx; j++)
        temp1[j] = cabsl(psi[j]);

    arr_mult_l(temp1, temp1, temp1, nx);   // |psi|^2
    arr_mult_l(temp2, temp1, temp1, nx);   // |psi|^4
    pe[m] = -0.5*arr_sum_l(temp2, nx)/arr_sum_l(temp1, nx);   // U

    for (int j = 0; j < nx; j++)
    {
        psi_r[ind2(m,j)] = creall(psi[j]);
        psi_i[ind2(m,j)] = cimagl(psi[j]);
    }	

    // Calculate KE
    fftwl_execute(f0_plan);                             // Get step's spectrum 
    for (int j = 0; j < nx; j++)
        spectrum[ind2(m,j)] = cabsl(psi_f[j])/nx;   // Save spectrum

    arr_mult_l(temp1, &spectrum[ind2(m, 0)], &spectrum[ind2(m,0)], nx);             //|psi_f|^2
    arr_mult_l(temp2, temp1, k2, nx);                //|psi_f|^2*k^2
    ke[m] =  0.5*arr_sum_l(temp2, nx)/arr_sum_l(temp1, nx);  // T   

    // Calculate Energy
    E[m] = pe[m] + ke[m];
    dE[m] = E[m] - E[0];
}
