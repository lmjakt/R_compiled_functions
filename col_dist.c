#include <R.h>
#include <Rinternals.h>
#include <math.h>

// a function that computes the euclidean distance between columns and
// returns a mirrored matrix. I have a preference for this over the
// behaviour of the 'dist' function in R as I can then dump the resulting
// matrix directly into 'image' to visualise it. 

// note that this function does not do any normalising or otherwise take into
// account the lengths of the vectors.. 

double distance( double *a, double *b, int l ){
  double dist = 0;
  for(int i = 0; i < l; ++i )
    dist += (a[i] - b[i]) * (a[i] - b[i]);
  return( sqrt(dist) );
}

// this accepts only a matrix of doubles
SEXP col_dist( SEXP r_m ){
  if( TYPEOF(r_m) != REALSXP )
    error("col_dist: data must be a matrix of real values\n");
  SEXP r_dims = PROTECT( getAttrib( r_m, R_DimSymbol ));
  if(length( r_dims ) != 2 )
    error("col_dist: data must be a matrix ndim %d\n", length( r_dims ));
  int *dims = INTEGER( r_dims );
  double *m = REAL( r_m );
  int nrow = dims[0];
  int ncol = dims[1];
  UNPROTECT(1);

  // we then set up the returning matrix
  SEXP r_dist = PROTECT( allocMatrix( REALSXP, ncol, ncol ) );
  double *dist = REAL( r_dist );

  for(int i=0; i < ncol; ++i){
    // set the diagonal; we do not need to calculate this as it 
    // should always be 0
    dist[ i * ncol + i ] = 0;
    for(int j=i+1; j < ncol; ++j){
      double d = distance( m + i * nrow, m + j * nrow, nrow );
      // this needs to be put into position i,j and j,i
      dist[ i * ncol + j ] = d;
      dist[ j * ncol + i ] = d;
    }
  }
  // which should be everything.
  UNPROTECT(1);
  return( r_dist );

}

