#ifndef SCRAMBLED_SOBOL_H_
#define SCRAMBLED_SOBOL_H_

/* The following max dim definition is for the top mass library only.
 * This is not a limitation of the Sobol algorithm implementation.
 */

#ifdef __cplusplus
extern "C" {
#endif

void i4_scrambled_sobol ( int dim_num, int *seed, float quasi[ ] );
void i8_scrambled_sobol ( int dim_num, long int *seed, double quasi[ ] );
  
#ifdef __cplusplus
}
#endif

#endif /* SCRAMBLED_SOBOL_H_ */
