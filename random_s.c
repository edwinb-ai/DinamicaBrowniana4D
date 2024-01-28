#include <math.h>
#include "random_s.h"

void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        mt[mti] &= 0xffffffffUL;
    }
}

double genrand_real1(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generation de N mots a un temp */
        int kk;

        if (mti == N+1)   /* si init_genrand() n'a pas êté apellé, */
            init_genrand(5489UL); /* un grain initiales par depit est utilisé */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y*(1.0/4294967296.0); 
    /* division pour 2^32-1 */ 
}

double gasdev(void)
{
    static double v_1;
    static double v_2;
    double rsq;
    static double g_dev;
    static int iset = 0;
    static double gset;

    if (iset == 0)
    {
        v_1 = 2.0*genrand_real1()-1.0;
        v_2 = 2.0*genrand_real1()-1.0;
        rsq = v_1*v_1 + v_2*v_2;

        while ((rsq >= 1.0) || (rsq == 0.0))
        {
            v_1 = 2.0*genrand_real1()-1.0;
            v_2 = 2.0*genrand_real1()-1.0;
            rsq = v_1*v_1 + v_2*v_2;
        }
        rsq = sqrt(-2.0*log(rsq)/rsq);
        gset = v_1*rsq;
        g_dev = v_2*rsq;
        iset = 1;
    }
    else
    {
        g_dev = gset;
        iset = 0;
    }
    return g_dev;
}