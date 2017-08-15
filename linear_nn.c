#include <R.h>
#include <Rinternals.h>

// given a matrix of all against all self distances, creates a 
// nearest neighbourhood network of neighbours within.
// r_delta should a sqaure matrix of all against all distances

SEXP linear_nn(SEXP r_delta, SEXP r_maxSkip)
{
  if( TYPEOF(r_delta) != REALSXP )
    error("The distances must be given as real values\n");
  if( TYPEOF(r_maxSkip) != INTSXP )
    error("The max skip variable should be an integer\n");
  SEXP r_dims = PROTECT( getAttrib( r_delta, R_DimSymbol ));
  if(length(r_dims) != 2)
    error("The distances must be expresses as a 2 dimensional (square) matrix");
  int *dims = (int*)INTEGER( r_dims );
  if( dims[0] != dims[1] )
    error("The distance matrix must be square not: %d x %d\n", dims[0], dims[1]);
  if( dims[0] < 2 )
    error("Less than two entries in the distance matrix. This is not what you wanted\n");

  ssize_t maxSkip = asInteger( r_maxSkip );
  if(maxSkip < 0)
    error("Negative maxSkip does not make sense\n");
  
  double *delta = REAL(r_delta);

  // The questions then is how to express sets of neighbours... 
  // the simplest would seem to be to simply define the relationships
  // using pointers and then go through and assign cluster ids.
  // That way we do not need to resize data. 
  
  // keep an index rather than ponters since we do not have any
  // complicated data structure, it makes more sense to simply
  // have a single link.
  size_t *links = malloc( sizeof(size_t*) * dims[0] );
  
  for(ssize_t i=0; i < (ssize_t)dims[0]; ++i){
    // find the nearest neighbour by looking at the distance up and dowstream
    // of itself..
    ssize_t b = i - maxSkip < 0 ? 0 : i - maxSkip;
    ssize_t e = i + maxSkip < dims[0] ? i + maxSkip : dims[0] - 1;
    ssize_t offset = i * dims[0];
    links[i] = 0; // default value
    ssize_t min_j = (b == i) ? e : b;
    ssize_t min_dist = (b == i) ? delta[ offset + e ] : delta[ offset + b ];
    for(ssize_t j=b; j <= e; ++j){
      if( j == i ) continue;
      if(min_dist > delta[ offset + j ]){
	min_dist = delta[ offset + j ];
	min_j = j;
      }
    }
    links[ i ] = (size_t)min_j;
  }
  
  // To harvest the groups we need to go through and follow the links until they are completely
  // self-contained; i.e. points to a previously assigned member. Not really sure how well this
  // will work
  
  // count membership from 1; that way we can leave 0 for non-assigned entries.. 
  SEXP r_membership = PROTECT( allocVector( INTSXP, dims[0] ));
  int *membership = INTEGER( r_membership );
  bzero( (void*)membership, sizeof( int ) * dims[0] );

  int current_cluster = 1;
  for(size_t i=0; i < dims[0]; ++i){
    if(membership[i]){
      // assign membership to closest neighbour if not already assigned
      if( !membership[ links[i] ] )
	membership[ links[i] ] = membership[i];
      continue;
    }

    if( membership[ links[i] ] ){
      membership[ i ] = membership[ links[i] ];
    }else{
      membership[ i ] = current_cluster;
      membership[ links[i] ] = current_cluster;
      ++current_cluster;
    }
  }
  free( links );
  UNPROTECT(2);
  return( r_membership );
}
