#include <math.h>
#include "random_s.h"

double norm(void) {
    static double v_1;
    static double v_2;
    double rsq;
    static double g_dev;
    static int iset = 0;
    static double gset;

    if (iset == 0) {
        v_1 = 2.0*genrand_real2()-1.0;
        v_2 = 2.0*genrand_real2()-1.0;
        rsq = v_1*v_1 + v_2*v_2;

        while ((rsq >= 1.0) || (rsq == 0.0)) {
            v_1 = 2.0*genrand_real2()-1.0;
            v_2 = 2.0*genrand_real2()-1.0;
            rsq = v_1*v_1 + v_2*v_2;
        }
        rsq = sqrt(-2.0*log(rsq)/rsq);
        gset = v_1*rsq;
        g_dev = v_2*rsq;
        iset = 1;
    }
    else {
        g_dev = gset;
        iset = 0;
    }
    
    return g_dev;
}