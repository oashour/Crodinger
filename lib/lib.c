#include "lib.h"

void simulation_destroy(struct simulation *sim)
{
    assert(sim != NULL);

    free(sim);
}

void simulation_print(struct simulation *sim)
{
    printf("Grid parameters: \n");
    printf("Temporal spacing     (dt): %.4"FS"\n", sim->dt);
    printf("Spatial spacing      (dx): %.4"FS"\n", sim->dx);
    printf("# x nodes            (nx): %d  \n",  sim->nx);
    printf("# t nodes            (nt): %d  \n",  sim->nt);
    printf("Max time             (tm): %.2"FS"\n", sim->tm);
    printf("Box size             (l) : %.4"FS"\n", sim->l);
    printf("\n");

    printf("Algorithmic information: \n");
    printf("Algorithm:           (--): %d%c\n", sim->order, sim->type);
    printf("\n");
        
    printf("Initial wave function: \n");
    printf("Modulation Amplitude (A1): %.4"FS"\n", sim->A1);
    printf("Constant Background  (A0): %.4"FS"\n", sim->A0);
    printf("Fundamental Freq     (Om): %.4"FS"\n", sim->Omega);
    printf("AB Parameter         (a) : %.4"FS"\n", sim->q);
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

void print_output(MYTYPE *psi_i, MYTYPE *psi_r, MYTYPE *x, 
                  MYTYPE *spectrum, MYTYPE *ke, MYTYPE *pe, 
                  MYTYPE *E, MYTYPE *dE, struct simulation *sim)
{
    FILE *fp = fopen(sim->f_param, "w");

    fprintf(fp, "%.13"FS"\n", sim->dt);
    fprintf(fp, "%d\n",    sim->nx);
    fprintf(fp, "%.13"FS"\n", sim->tm);
    fprintf(fp, "%.13"FS"\n", sim->A1);
    fprintf(fp, "%.13"FS"\n", sim->q);
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
    MYTYPE dt, tm, l; 
    MYTYPE A1 = 0, q = 0;
    int order, nx, sampling, l_mult, initial; 
    char an_choice;
    char type;
    MYTYPE nu = 0;
    char f_psi_i[256], f_psi_r[256], f_spectrum[256], f_param[256], f_x[256];
    char f_pe[256], f_ke[256], f_E[256], f_dE[256];
    
    // Basic grid parameters
    printf("dt: ");                         // Grid temporal spacing          
    scanf("%"FS, &dt);
    printf("nx: ");                         // Number of Fourier modes
    scanf("%d", &nx);                       
    printf("tm: ");                         // Maximum time
    scanf("%"FS, &tm);

    // Initial Conditions
    printf("For the initial wavefunction, [1] denotes "
           "Peregrine soliton, [2] denotes background.\n"
           "Initial wavefunction code: ");
    scanf("%d", &initial);
    if (initial == BACKGROUND)
    {
        printf("A1: ");                         // Cosine amplitude 
        scanf("%"FS, &A1);
        printf("Do you want to pick a or n? [a/n]");                          // Breather parameter
        scanf(" %c", &an_choice);
        if (an_choice == 'a')
        {
            printf("a: ");                          // Breather parameter
            scanf("%"FS, &q);
        }
        else if (an_choice == 'n')
        {
            printf("n: ");                          // Breather parameter
            scanf("%"FS, &nu);
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
        scanf("%"FS, &l);
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
    const MYTYPE dx = (l/nx);		   // Spatial step size
    MYTYPE Omega = 0, A0 = 0;
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

void arr_mult_c(MYFFTW(_complex) *out, MYFFTW(_complex) *in1, MYFFTW(_complex) *in2, int size)
{
    for(int i =0; i < size; i++)
        out[i] = in1[i]*in2[i];
}

void arr_mult_l(MYTYPE *out, MYTYPE *in1, MYTYPE *in2, int size)
{
    for(int i =0; i < size; i++)
        out[i] = in1[i]*in2[i];
}

void arr_mult(MYFFTW(_complex) *psi, MYTYPE mult, int size)
{
    for (int i = 0; i < size; i++)
        psi[i] = psi[i]*mult;
}

MYTYPE arr_sum_l(MYTYPE *in, int size)
{
    MYTYPE out = 0;
    for (int i =0; i < size; i++)
        out += in[i];

    return out;
}

void arr_add_l(MYTYPE *out, MYTYPE *in1, MYTYPE *in2, int size)
{
    for (int i = 0; i < size; i++)
        out[i] = in1[i]+in2[i];
}

void arr_add(MYFFTW(_complex) *out, MYFFTW(_complex) *in1, MYFFTW(_complex) *in2, int size)
{
    for (int i = 0; i < size; i++)
        out[i] = in1[i]+in2[i];
}

void T2(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi_o, int nt, MYTYPE dt, MYTYPE *k2, int nx, MYFFTW(_plan) forward, MYFFTW(_plan) backward) 
{
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabsl(psi[i]) * cabsl(psi[i]) * dt/D(2.0))*psi[i];
    // Forward transform
    MYFFTW(_execute_dft)(forward, psi_o, psi_o);
    // Solve linear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(-I * k2[i]/D(2.0) * dt)*psi_o[i];
    // Backward transform
    MYFFTW(_execute_dft)(backward, psi_o, psi_o);
    // Normalize the transform
    arr_mult(psi_o, D(1.0)/nx, nx);
    // Solve nonlinear part
    for(int i = 0; i < nx; i++)
        psi_o[i] = cexp(I * cabsl(psi_o[i]) * cabsl(psi_o[i]) * dt/D(2.0))*psi_o[i];
}

void T4_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, int nt, MYTYPE dt, 
          MYTYPE *k2, int nx, MYFFTW(_plan) forward, MYFFTW(_plan) backward)
{
	    T2(psi, psi1, nt, dt/D(2.0), k2, nx, forward, backward);
	    T2(psi1, psi1, nt, dt/D(2.0), k2, nx, forward, backward);
        arr_mult(psi1, D(4.0)/D(3.0), nx);

	    T2(psi, psi2, nt, dt, k2, nx, forward, backward);
        arr_mult(psi2, -D(1.0)/D(3.0), nx);

        arr_add(psi, psi1, psi2, nx);
}

void T6_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, MYFFTW(_complex) *psi3,
          int nt, MYTYPE dt, MYTYPE *k2, int nx, MYFFTW(_plan) forward, MYFFTW(_plan) backward)
{
        T2(psi, psi1, nt, dt/D(3.0), k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/D(3.0), k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/D(3.0), k2, nx, forward, backward);
        arr_mult(psi1, D(81.0)/D(40.0), nx);

        T2(psi, psi2, nt, dt/D(2.0), k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/D(2.0), k2, nx, forward, backward);
        arr_mult(psi2, D(-16.0)/D(15.0), nx);

        T2(psi, psi3, nt, dt, k2, nx, forward, backward);
        arr_mult(psi3, D(1.0)/D(24.0), nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
}

void T8_M(MYFFTW(_complex) *psi, MYFFTW(_complex) *psi1, MYFFTW(_complex) *psi2, MYFFTW(_complex) *psi3,
          MYFFTW(_complex) *psi4, int nt, MYTYPE dt, MYTYPE *k2, int nx, MYFFTW(_plan) forward, 
          MYFFTW(_plan) backward)
{
        T2(psi, psi1, nt, dt/D(4.0), k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/D(4.0), k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/D(4.0), k2, nx, forward, backward);
        T2(psi1, psi1, nt, dt/D(4.0), k2, nx, forward, backward);
        arr_mult(psi1, D(1024.0)/D(315.0), nx);

        T2(psi, psi2, nt, dt/D(3.0), k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/D(3.0), k2, nx, forward, backward);
        T2(psi2, psi2, nt, dt/D(3.0), k2, nx, forward, backward);
        arr_mult(psi2, D(-729.0)/D(280.0), nx);

        T2(psi, psi3, nt, dt/D(2.0), k2, nx, forward, backward);
        T2(psi3, psi3, nt, dt/D(2.0), k2, nx, forward, backward);
        arr_mult(psi3, D(16.0)/D(45.0), nx);

        T2(psi, psi4, nt, dt, k2, nx, forward, backward);
        arr_mult(psi4, D(-1.0)/D(360.0), nx);

        arr_add(psi, psi1, psi2, nx);
        arr_add(psi, psi, psi3, nx);
        arr_add(psi, psi, psi4, nx);
}

void T4_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx, 
          MYFFTW(_plan) forward, MYFFTW(_plan) backward) 
{
        MYTYPE s = pow(2, D(1.0)/D(3.0));
        MYTYPE os = D(1.0)/(D(2.0)-s);

        MYTYPE ft = os;
        MYTYPE bt = -s*os;

        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, bt*dt, k2, nx, forward, backward);
        T2(psi, psi, nt, ft*dt, k2, nx, forward, backward);
}

void T6_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx, 
          MYFFTW(_plan) forward, MYFFTW(_plan) backward) 
{
        MYTYPE s = pow(2, D(1.0)/D(5.0));
        MYTYPE os = D(1.0)/(D(2.0)-s);

        MYTYPE ft = os;
        MYTYPE bt = -s*os;

        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T4_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T4_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

void T8_S(MYFFTW(_complex) *psi, int nt, MYTYPE dt, MYTYPE *k2, int nx, 
          MYFFTW(_plan) forward, MYFFTW(_plan) backward) 
{
        MYTYPE s = pow(2, D(1.0)/D(7.0));
        MYTYPE os = D(1.0)/(D(2.0)-s);

        MYTYPE ft = os;
        MYTYPE bt = -s*os;

        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
        T6_S(psi, nt, bt*dt, k2, nx, forward, backward);
        T6_S(psi, nt, ft*dt, k2, nx, forward, backward);
}

void energy(MYFFTW(_complex) *psi, MYTYPE *psi_r, MYTYPE *psi_i, MYTYPE *spectrum,
            MYTYPE *ke, MYTYPE *pe, MYTYPE *E, MYTYPE *dE,
            MYFFTW(_plan) f0_plan, MYFFTW(_complex) *psi_f, MYTYPE *k2, int m, int nx)
{
    MYTYPE *temp1 = (MYTYPE*)malloc(sizeof(MYTYPE) * nx);
    MYTYPE *temp2 = (MYTYPE*)malloc(sizeof(MYTYPE) * nx);

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
    MYFFTW(_execute)(f0_plan);                             // Get step's spectrum 
    for (int j = 0; j < nx; j++)
        spectrum[ind2(m,j)] = cabsl(psi_f[j])/nx;   // Save spectrum

    arr_mult_l(temp1, &spectrum[ind2(m, 0)], &spectrum[ind2(m,0)], nx);             //|psi_f|^2
    arr_mult_l(temp2, temp1, k2, nx);                //|psi_f|^2*k^2
    ke[m] =  0.5*arr_sum_l(temp2, nx)/arr_sum_l(temp1, nx);  // T   

    // Calculate Energy
    E[m] = pe[m] + ke[m];
    dE[m] = abs(E[m] - E[0]);
}
