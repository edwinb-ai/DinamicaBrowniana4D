#ifndef FUNCTIONSL_H
#define FUNCTIONSL_H

// Algunas variables globales
static const double lambda = 50.0;
static const double a_param = 134.5526623421209;
static const double b_param = 1.0204081632653061;
static const double temp = 1.548382;
static const double mt_n = 100000;
static const int mp = 2500;
static const int nm = 2500;
// static const double pisq = 3.141592653589793;

// Funciones generales del código
void iniconf(double* x, double* y, double* z, double* w, const double rho,
const double rc, const int num_part);

double rdf_force(double* x, double* y, double* z, double*w, double* fx, double* fy, double* fz,
double* fw, const int num_part, const double box_l);

double hardsphere(double rij);

void position(double* x, double* y, double* z, double* w, double* fx, double* fy, double* fz, double* fw, const double dtt,
const double box_l, const int num_part, const int pbc);

void gr(double* x, double* y, double* z, double* w, double* g, const int num_part, const double box_l);

void difusion(const int nprom, const int n_part, double* cfx, double* cfy, double* cfz,
double* cfw, double* wt);

#endif