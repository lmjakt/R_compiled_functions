#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <math.h>

// given a set of values and a set of mid points, calculates a kernel
// smoothened histogram for the mid points using the specified kernel
// and the distance scale for the kernel (as the distance of one unit of the kernel)

// note that the kernel should be symmetric and only half given such that the distance
// between the mid-point and value can be used to determine the kernel multipland.
// Note that the function does nothing to compensate for edge effects as it handles
// sparse data and it is not clear how to do this most effectively... 
SEXP smoothed_hist(SEXP r_values, SEXP r_mids, SEXP r_kernel, SEXP r_k_unit){
  if( !isReal( r_values ) || !isReal( r_mids ) || !isReal( r_kernel ) || !isReal(r_k_unit) )
    error("All arguments should be real vectors\n");

  double k_unit = asReal( r_k_unit );
  double *values = REAL(r_values);
  double *mids = REAL(r_mids);
  double *kernel = REAL( r_kernel );

  size_t values_l = length( r_values );
  size_t mids_l = length( r_mids );
  size_t kernel_l = length( r_kernel );
  
  if( values_l < 1 || mids_l < 1 || kernel_l < 1 )
    error("Empty values, mids or kernel not allowed\n");
  
  if( k_unit <= 0 )
    error("Positive k_unit required\n");
 // lets allocate some space..
  SEXP r_counts = PROTECT(allocVector( REALSXP, mids_l ));
  double *counts = REAL( r_counts );
  bzero( (void*)counts, sizeof(double) * mids_l );

  // this is a bit ugly, if the values were sorted then one could do something
  // with them, but never mind..
  size_t o = 0;

  for(size_t i=0; i < values_l; ++i){
    for(size_t j=0; j < mids_l; ++j){
      o = (size_t) ( fabs(values[i] - mids[j]) / k_unit );
      if(o < kernel_l)
	counts[j] += kernel[ o ];
    }
  }
  UNPROTECT(1);
  return( r_counts );
}
