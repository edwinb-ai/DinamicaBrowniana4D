#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "random_s.h"
#include "functionsl.h"

void iniconf(double* x, double* y, double* z, double* w, double rho, double rc, int num_part)
{
    // Definir la distancia segun la densidad
    double dist = pow(rho, -1.0/4.0);

    // Inicializar las primeras posiciones
    x[0] = -rc + (dist / 2.0);
    y[0] = -rc + (dist / 2.0);
    z[0] = -rc + (dist / 2.0);
    w[0] = -rc + (dist / 2.0);

    for (int i = 0; i < num_part; i++)
    {
        x[i+1] = x[i] + dist;
        y[i+1] = y[i];
        z[i+1] = z[i];
        w[i+1] = w[i];

        if (x[i+1] > rc)
        {
            x[i+1] = -rc + (dist / 2.0);
            y[i+1] += dist;
            z[i+1] = z[i];
            w[i+1] = w[i];

            if (y[i+1] > rc)
            {
                x[i+1] = -rc + (dist / 2.0);
                y[i+1] = -rc + (dist / 2.0);
                z[i+1] += dist;
                w[i+1] = w[i];

                if (z[i+1] > rc)
                {
                    x[i+1] = -rc + (dist / 2.0);
                    y[i+1] = -rc + (dist / 2.0);
                    z[i+1] = -rc + (dist / 2.0);
                    w[i+1] += dist;
                }
            }
        }
    }
}

double hardsphere(double r_pos)
{
    double uij = 0.0;

    uij = (a_param/temp) * (pow(1.0/r_pos, lambda) - pow(1.0/r_pos, lambda-1.0));

    uij += 1.0 / temp;

    return uij;
}

double rdf_force(double* x, double* y, double* z, double*w, double* fx, double* fy, double* fz,
double* fw, const int num_part, const double box_l)
{
    // Parámetros
    double rc = box_l/2.0;
    double d_r = rc / nm;

    // Inicializar algunas variables de la posicion
    double xij = 0.0, yij = 0.0, zij = 0.0, rij = 0.0;
    double fij = 0.0, wij = 0.0;
    double uij = 0.0, ener = 0.0;
    int nbin = 0;
    size_t i = 0, j = 0;

    // Inicializar arreglos para la fuerza
    for (i = 0; i < num_part; i++)
    {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
        fw[i] = 0.0;
    }

    // #pragma omp parallel for num_threads(3) default(shared) private(xij,yij,zij,wij,i,j,rij) reduction(+:ener)
    for (i = 0; i < num_part; i++)
    {
        for (j = i+1; j < num_part-1; j++)
        {
            // Siempre inicializar en cero
            uij = 0.0;  
            fij = 0.0;  

            // Contribucion de pares
            xij = x[i] - x[j];
            yij = y[i] - y[j];
            zij = z[i] - z[j];
            wij = w[i] - w[j];

            // Condiciones de frontera
            xij -= (box_l * round(xij/box_l));
            yij -= (box_l * round(yij/box_l));
            zij -= (box_l * round(zij/box_l));
            wij -= (box_l * round(wij/box_l));

            rij = sqrt(xij*xij + yij*yij + zij*zij + wij*wij);
            if (rij < rc)
            {
                // Siempre se calcula la fuerza
                if (rij < b_param)
                {
                    uij = hardsphere(rij);
                    fij = lambda*pow(1.0/rij, lambda+1.0) - (lambda-1.0)*pow(1.0/rij, lambda);
                    fij *= (a_param/temp);
                }
                else
                {
                    uij = 0.0;
                    fij = 0.0;
                }
                // Actualizar los valores de las fuerzas
                fx[i] += (fij*xij)/rij;
                fy[i] += (fij*yij)/rij;
                fz[i] += (fij*zij)/rij;
                fw[i] += (fij*wij)/rij;

                fx[j] -= (fij*xij)/rij;
                fy[j] -= (fij*yij)/rij;
                fz[j] -= (fij*zij)/rij;
                fw[j] -= (fij*wij)/rij;
                ener = ener + uij;
            }
        }
    }

    return ener;
}

void position(double* x, double* y, double* z, double* w, double* fx, double* fy,
double* fz, double* fw, const double dtt,
const double box_l, const int num_part, const int pbc)
{
    // Inicializar algunas variables
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double dw = 0.0;
    double sigma = sqrt(2.0*dtt);

    for (int i = 0; i < num_part; i++)
    {
        dx = sigma * norm();
        dy = sigma * norm();
        dz = sigma * norm();
        dw = sigma * norm();

        x[i] += fx[i]*dtt + dx;
        y[i] += fy[i]*dtt + dy;
        z[i] += fz[i]*dtt + dz;
        w[i] += fw[i]*dtt + dw;

        if (pbc == 1)
        {
            x[i] -= (box_l * round(x[i]/box_l));
            y[i] -= (box_l * round(y[i]/box_l));
            z[i] -= (box_l * round(z[i]/box_l));
            w[i] -= (box_l * round(w[i]/box_l));
        }
    }
}

void gr(double* x, double* y, double* z, double* w, double* g,
const int num_part, const double box_l)
{
    // Parámetros
    double rc = box_l/2.0;
    double d_r = rc / nm;

    int nbin = 0;
    int i = 0, j = 0;
    double xij = 0.0, yij = 0.0, zij = 0.0, rij = 0.0;
    double wij = 0.0;

    for (i = 0; i < num_part; i++)
    {
        for (j = i+1; j < num_part-1; j++)
        {

            // Contribucion de pares
            xij = x[j] - x[i];
            yij = y[j] - y[i];
            zij = z[j] - z[i];
            wij = w[j] - w[i];

            // Condiciones de frontera
            xij -= (box_l * round(xij/box_l));
            yij -= (box_l * round(yij/box_l));
            zij -= (box_l * round(zij/box_l));
            wij -= (box_l * round(wij/box_l));

            rij = sqrt(xij*xij + yij*yij + zij*zij + wij*wij);

            if (rij < rc)
            {
                nbin = (int)(rij/d_r) + 1;
                if (nbin <= nm)
                {
                    g[nbin] += 2.0;
                }
            }
        }
    }
}

void difusion(const int nprom, const int n_part, double* cfx, double* cfy, double* cfz,
double* cfw, double* wt)
{
    double dif2 = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0, dw = 0.0;
    double aux2 = 0.0;

    for (int i = 0; i < nprom; i++)
    {
        dif2 = 0.0;
        for (int j = 0; j < nprom-i; j++)
        {
            for (int k = 0; k < n_part; k++)
            {
                dx = cfx[(j+i)*mp + k] - cfx[j*mp + k];
                dy = cfy[(j+i)*mp + k] - cfy[j*mp + k];
                dz = cfz[(j+i)*mp + k] - cfz[j*mp + k];
                dw = cfw[(j+i)*mp + k] - cfw[j*mp + k];
                dif2 += dx*dx + dy*dy + dz*dz + dw*dw;
            }
        }
        aux2 = n_part * (nprom - i);
        wt[i] += (dif2/aux2);
    }
}