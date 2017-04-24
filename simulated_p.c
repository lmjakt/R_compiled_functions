#include <R.h>
#include <Rinternals.h>

// this function takes two arguments:
// a set of test statistics from an experiment
// a set of test statistics obtained from a simulated
// situation where the null hypothesis is true.

// both vectors need to be sorted from small to large
// the p-values returned are the probability of observing
// a statistic larger than to the given one. Note that it's minimal
// value will be 1 / permutation number.

// The function returns a matrix containing p-values and
// the false discovery rates of those p-value.

// WARNING r_test and r_ctl must be sorted vectors..

SEXP simulated_p(SEXP r_test, SEXP r_ctl){
  if(!isVector(r_test) || !isVector(r_ctl))
    error("Both arguments should be vectors\n");
  if(!isReal(r_test) || !isReal(r_ctl))
    error("Both arguments should be real vectors\n");
  // it seems that we do have an isUnsorted function here
  // I am using 0 here to indicat the value of a Rboolean type. I should find
  // where it is defined and then see if there is a better way of defining it.
  if(isUnsorted(r_test, FALSE) || isUnsorted(r_ctl, FALSE))
    error("Both test and ctl values have to be sorted, low to high\n");

  int t_length = length(r_test);
  int c_length = length(r_ctl);
  double repNo = (double)c_length / (double)t_length;
  
  double *test = REAL(r_test);
  double *ctl = REAL(r_ctl);
  
  // allocate a matrix for the resulting p-values..
  // one column for p one column for fdr
  SEXP r_p = PROTECT( allocMatrix( REALSXP, t_length, 2) );
  double *p = REAL(r_p);
  double *fdr = REAL(r_p) + t_length;

  // this is a bit hacky, but should set the full set of values.. 
  memset(p, 0, t_length * 2 * sizeof(double));

  size_t c_i = 0;
  for(size_t t_i=0; t_i < t_length; ++t_i){
    while( c_i < c_length && ctl[ c_i ] <= test[ t_i ])
      ++c_i;
    // the number of values which are larger than or equal to the given value
    // under the null hypothesis (or the full data set including the test statistics)
    p[t_i] = (double)(t_length + c_length - (c_i + t_i)) / (double)(t_length + c_length);
    
    // the false discovery rate; the ratio of the rate of discovery under
    // the null hypothesis to the rate of discovery under the test situation
    fdr[t_i] = (double)(c_length - c_i) / ((double)(t_length - t_i) * repNo);
  }
  UNPROTECT(1);
  return( r_p );
}

