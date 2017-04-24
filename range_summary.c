#include <R.h>
#include <Rinternals.h>

// Given a range defined by
// a matrix of starts and ends
// and a set of scores return summary scores for each range

// use reals for everything so that we can handle more than 2e9 offsets


SEXP range_summary(SEXP r_start, SEXP r_end, SEXP r_scores){
  if(!isVector(r_scores) || !isVector(r_start) || !isVector(r_end))
    error("all arguments should be vectors\n");
  if(!isReal(r_start) || !isReal(r_end) || !isReal(r_scores))
    error("all arguments should be Real values\n");

  size_t rN = length(r_start);
  if(length( r_end ) != rN || rN  == 0)
    error("Length of starts and ends must be the same\n");
  
  // calculate the following for each region:
  // 1. mean
  // 2. min
  // 3. max
  // 4. width
  // 5. sum  (not really necessary, but we can put it there in any case)
  
  // allocate the results structure
  // NOTE that it would probably be faster to make a transposed
  // matrix, but it's unlikely to be much of an issue (and this looks easier)
  SEXP r_summaries = PROTECT(allocMatrix( REALSXP, rN, 5 ));
  double *summary = REAL(r_summaries);
  double *mean = summary;
  double *min = summary + rN;
  double *max = summary + 2 * rN;
  double *width = summary + 3 * rN;
  double *sum = summary + 4 * rN;

  size_t sN = length( r_scores );
  double *scores = REAL( r_scores );

  double *start = REAL( r_start );
  double *end = REAL( r_end );

  for(size_t i=0; i < rN; ++i){
    size_t b = (size_t)start[i];
    size_t e = (size_t)end[i];
    if( e < b || e >= sN ){
      warning("range out of bounds or negative\n");
      mean[i] = 0;
      min[i] = 0;
      max[i] = 0;
      width[i] = 0;
      sum[i] = 0;
      continue;
    }
    sum[i] = scores[ b ];
    min[i] = max[i] = scores[ b ];
    size_t j = b + 1;
    while(j <= e){
      sum[i] += scores[j];
      max[i] = max[i] < scores[j] ? scores[j] : max[i];
      min[i] = min[i] > scores[j] ? scores[j] : min[i];
      ++j;
    }
    width[i] = 1 + e - b;
    mean[i] = sum[i] / width[i];
  }

  UNPROTECT(1);
  return( r_summaries );
}
