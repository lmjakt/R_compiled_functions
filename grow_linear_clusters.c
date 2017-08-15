#include <R.h>
#include <Rinternals.h>

// given a matrix of all against all distances grows square clusters until
// the maximum internal inter-node distance exceeds a threshold.
// Iniitally this will just start at the first one and go forwards. It seems
// that this might have order based effects, but I can't think of a better way
// to start with.

SEXP grow_linear_clusters( SEXP r_delta, SEXP r_max_delta )
{
  if( TYPEOF(r_delta) != REALSXP )
    error("The distances must be given as real values\n");
  if( TYPEOF(r_max_delta) != REALSXP )
    error("The maximum distance must be given as a real value\n");
  
  SEXP r_dims = PROTECT( getAttrib( r_delta, R_DimSymbol ));
  if( length(r_dims) != 2 )
    error("The distances must be given as a two-dimensional matrix\n");
  
  int *dims = INTEGER(r_dims);
  if(dims[0] != dims[1] || dims[0] == 0)
    error("The distance matrix must be square! Not: %d x %d\n", dims[0], dims[1]);

  int nrow=dims[0];

  double max_delta = asReal( r_max_delta );
  double *delta = REAL( r_delta );

  // represent cluster membership simply as a number for the given window. We can
  // consider later on to make a more complicated data structure containing more
  // information.

  SEXP r_clusters = PROTECT( allocVector( INTSXP, nrow ));
  int *clusters = INTEGER( r_clusters );
  bzero( (void*)clusters, sizeof(int) * nrow );

  // the begin and end of the current cluster
  size_t b, e, ee;
  // the maximum distance

  // we don't need to know the maximum delta, but in future versions we may wish
  // to return it as part of a list as this is useful information for the end user.
  // but first 
  // double max_d = 0; 
  int current_cluster = 1;

  for(b = 0; b < nrow; ++b){
    //    max_d = 0;
    for(e = b + 1; b < nrow; ++e){
      clusters[e] = current_cluster;
      for(ee = e-1; ee >= b; --ee){
	if( delta[ e * nrow + ee ] > max_delta ){
	  e = e - 1;
	  ee = -1;
	  clusters[e] = 0;
	  break;
	}
      }
      if(ee == -1)
	break;
    }
    ++current_cluster;
    b = e + 1;
  }
  UNPROTECT(2);
  return( r_clusters );
}
