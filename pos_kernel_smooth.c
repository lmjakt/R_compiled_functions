#include <R.h>
#include <Rinternals.h>

// Given a set of values at specified positions and a kernel, returns
// a kernel density estimate for each position.

// all arguments should be double vector

SEXP pos_kernel_smooth( SEXP r_values, SEXP r_positions, SEXP r_kernel ){
  if( !isVector(r_values) || !isVector(r_positions) || !isVector( r_kernel ))
    error("All arguments should be vectors\n");
  if( !isReal(r_values) || !isReal(r_positions) || !isReal(r_kernel) )
    error("All arguments should be real vectors\n");
  
  size_t vN = length(r_values);
  if( length( r_positions ) != vN )
    error("Positions and r_values must be of the same length\n");
  
  // finally, the kernel should be of an odd length and have a length of
  // at least 3; otherwise nonsense..
  size_t kN = length(r_kernel);
  if( !(kN % 2) || kN < 3 )
    error("The kernel must have a lenght that is odd and is larger than 3\n");
  
  // the positions must also be sorted for this to make any sense:
  if( isUnsorted( r_positions, FALSE ) )
    error("The positions must be sorted from low to high\n");
  
  // Then we can start
  double *values = REAL(r_values);
  double *positions = REAL(r_positions);
  double *kernel = REAL(r_kernel);

  SEXP r_dens = PROTECT(allocVector( REALSXP, vN ) );
  double *dens = REAL(r_dens);

  for(size_t d_i=0; d_i < vN; ++d_i){
    size_t beg = d_i;
    size_t end = d_i;

    // determine where we start on the kernel...
    while( beg > 1 && positions[ d_i ] - positions[ beg - 1 ] < (double)(kN / 2) )
      --beg;
    while( end < (vN - 2) && positions[ end + 1 ] - positions[ d_i ] < (double)(kN / 2) )
      ++end;

    dens[ d_i ] = 0;
    for(size_t v_i=beg; v_i <= end; ++v_i){
      size_t k_offset = (kN / 2) + (size_t)(positions[ v_i ] - positions[ d_i ]);
      dens[ d_i ] += kernel[ k_offset ] * values[ v_i ];
    }
  }
  
  UNPROTECT(1);
  return( r_dens );
}
