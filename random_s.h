#ifndef RANDOMS_H
#define RANDOMS_H

/* Parametres  du periode*/  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

void init_genrand(unsigned long s);
double genrand_real1(void);
double gasdev(void);

#endif