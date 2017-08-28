#ifndef SOBOL_H_
#define SOBOL_H_

/* The following max dim definition is for the top mass library only.
 * This is not a limitation of the Sobol algorithm implementation.
 */
#define MAX_SOBOL_DIM 40

#ifdef __cplusplus
extern "C" {
#endif

double d_uniform_01 ( int *seed );
int i4_bit_hi1 ( int n );
int i4_bit_lo0 ( int n );
void i4_sobol ( int dim_num, int *seed, float quasi[ ] );
int i4_uniform ( int b, int c, int *seed );
unsigned int i4_xor ( unsigned int i, unsigned int j );
int i8_bit_hi1 ( long int n );
int i8_bit_lo0 ( long int n );
void i8_sobol ( int dim_num, long int *seed, double quasi[ ] );
unsigned long int i8_xor ( unsigned long int i, unsigned long int j );
long int i8_uniform ( long int b, long int c, int *seed );

#ifdef __cplusplus
}
#endif

#endif /* SOBOL_H_ */
