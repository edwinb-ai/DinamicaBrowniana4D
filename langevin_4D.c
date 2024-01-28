#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "functionsl.h"
#include "pseudorng.h"

int main(int argc, char const *argv[]) {
    // Archivos para trabajar
    FILE *f_iniconf;
    FILE *f_gr;
    FILE *f_final;
    FILE *f_ener;
    FILE *f_wt;

    // pi al cuadrado
    const double pisq = M_PI*M_PI;

    // Numero de particulas
    int n_part = pow(4, 4);
    // Packing fraction
    double phi = atof(argv[1]);
    double rho = 32.0*phi/pisq;
    // Opciones de linea de comando
    int nct = atoi(argv[2]);
    unsigned int ncp = atoi(argv[3]);
    // Paso de tiempo
    double d_tiempo = atof(argv[4]);
    int thermal = atoi(argv[5]);

    // Tamaño de caja
    double l_caja = pow(n_part / rho, 1.0/4.0);
    double radio_c = l_caja/2.0;
    double dr = radio_c/nm;
    // double dq = M_PI/radio_c;

    // Mostrar información del sistema
    printf("El tamaño de la caja es: %f\n", l_caja);
    printf("Distance media entre partículas: %f\n", pow(rho, -1.0/4.0));
    printf("Radio de corte: %f\n", radio_c);

    // Inicializar el RNG, los siguientes comandos hacen que las semillas
    // siempre sean diferentes
    // Semilla para el PRNG
    unsigned long seed;
    // Abrimos el PRNG de Linux, tambien puede ser un numero fijo
    FILE *fp = fopen("/dev/urandom", "r");
    int tmp = fread(&seed, 1, sizeof(unsigned long), fp);
    if (tmp != sizeof(unsigned long)) printf ("error with seed\n");
    fclose(fp);
    // Mostramos la semilla para guardarla despues, por reproducibilidad
    printf("Semilla: %u\n", (int)(seed));
    init_genrand(seed);

    // Inicializar los arreglos
    double *x = calloc(mp, sizeof(double));
    double *y = calloc(mp, sizeof(double));
    double *z = calloc(mp, sizeof(double));
    double *w = calloc(mp, sizeof(double));
    double *fx = calloc(mp, sizeof(double));
    double *fy = calloc(mp, sizeof(double));
    double *fz = calloc(mp, sizeof(double));
    double *fw = calloc(mp, sizeof(double));
    double *g = calloc(nm, sizeof(double));
    double *t = calloc(mt_n, sizeof(double));
    double *h = calloc(nm, sizeof(double));
    double *wt = calloc(mt_n, sizeof(double));
    double *cfx = calloc(mt_n*mp, sizeof(double));
    double *cfy = calloc(mt_n*mp, sizeof(double));
    double *cfz = calloc(mt_n*mp, sizeof(double));
    double *cfw = calloc(mt_n*mp, sizeof(double));

    // Configuración inicial
    iniconf(x, y, z, w, rho, radio_c, n_part);
    f_iniconf = fopen("conf_inicial.dat", "w");
    for (int i = 0; i < n_part; i++) {
        fprintf(f_iniconf, "%.10f %.10f %.10f %.10f\n", x[i], y[i], z[i], w[i]);
    }
    fclose(f_iniconf);
    // Verificar que la energía es cero
    double ener = 0.0, enerpart;
    ener = rdf_force(x, y, z, w, fx, fy, fz, fw, n_part, l_caja);
    enerpart = ener / ((double)(n_part));
    printf("E/N: %.10f\n", enerpart);

    // Termalizar el sistema
    // Se puede utilizar una termalización anterior
    if (thermal == 1) {
        f_final = fopen("final_conf.dat", "r");
        for (int i = 0; i < n_part; i++) {
            fscanf(f_final, "%lf", &x[i]);
            fscanf(f_final, "%lf", &y[i]);
            fscanf(f_final, "%lf", &z[i]);
            fscanf(f_final, "%lf", &w[i]);
            fscanf(f_final, "%lf", &fx[i]);
            fscanf(f_final, "%lf", &fy[i]);
            fscanf(f_final, "%lf", &fz[i]);
            fscanf(f_final, "%lf", &fw[i]);
        }
        fclose(f_final);
    }
    // O crear una nueva
    else {
        f_ener = fopen("energia.dat", "w");
        f_final = fopen("final_conf.dat", "w");

        for (int i = 0; i < nct; i++) {
            position(x, y, z, w, fx, fy, fz, fw, d_tiempo, l_caja, n_part, 1);
            ener = rdf_force(x, y, z, w, fx, fy, fz, fw, n_part, l_caja);
            enerpart = ener / ((double)(n_part));
            if (i%10000 == 0) {
                printf("%d %.10f Thermal\n", i, enerpart);
            }
            if (i%100 == 0) {
                fprintf(f_ener, "%d %.10f\n", i, enerpart);
            }
        }
        fclose(f_ener);

        // Guardar la configuración final después de termalizar
        for (int i=0; i < n_part; i++) {
            fprintf(f_final, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
            x[i], y[i], z[i], w[i], fx[i], fy[i], fz[i], fw[i]);
        }
        fclose(f_final);   
    }

    // Calcular la g(r)
    int nprom = 0;
    for (int i = 0; i < ncp; i++) {
        position(x, y, z, w, fx, fy, fz, fw, d_tiempo, l_caja, n_part, 0);
        ener = rdf_force(x, y, z, w, fx, fy, fz, fw, n_part, l_caja);
        enerpart = ener / ((double)(n_part));
        if (i%10000 == 0) {
            printf("%d %.10f Average\n", i, enerpart);
        }
        if (i%10 == 0) {
            t[nprom] = d_tiempo*10.0*nprom;
            for (int j = 0; j < n_part; j++) {
                cfx[nprom*mp + j] = x[j];
                cfy[nprom*mp + j] = y[j];
                cfz[nprom*mp + j] = z[j];
                cfw[nprom*mp + j] = w[j];
            }
            nprom++;
            gr(x, y, z, w, g, n_part, l_caja);
        }
    }
    // Leer el nombre del archivo de la linea de comandos
    f_gr = fopen(argv[6], "w");
    double *r = calloc(nm, sizeof(double));
    double dv = 0.0;
    // double hraux = 0.0;
    for (int i = 1; i < nm; i++) {
        r[i] = (i-1)*dr;
        dv = dr*2.0*pisq*pow(r[i], 3);
        g[i] *= pow(l_caja, 4) / (pow(n_part, 2)*nprom*dv);
        h[i] = g[i] - 1.0;
        fprintf(f_gr, "%.10f %.10f %.10f\n", r[i], g[i], h[i]);
    }
    fclose(f_gr);

    // Llamar la función para el desplazamiento cuadrático medio
    difusion(nprom, n_part, cfx, cfy, cfz, cfw, wt);
    f_wt = fopen(argv[7], "w");
    for (int i = 0; i < ncp/10; i++) {
        fprintf(f_wt, "%.10f %.10f\n", t[i], wt[i]);
    }
    fclose(f_wt);

    // Liberar arreglos
    free(x); free(y); free(z); free(fx); free(fy); free(fz);
    free(w); free(fw), free(h), free(wt), free(t);
    free(cfx), free(cfy), free(cfz), free(cfw);
    free(r); free(g);

    return 0;

}