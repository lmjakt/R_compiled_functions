#include <R.h>
#include <Rinternals.h>

// This codes is inspired by cpg_distribute the code written by:
// Gos Micklem + Tim Cutts (2004)
//
// to define CpG islands in genomes.
// However rather than take a sequence of nucleotides it takes
// a sequence of scores and positions associated with those scores
// as well as a distance penalty...
// 
// Give the locations of CpGs as scores of CPGSCORE
// and a seperation penalty of 1 / position it should initially find the
// same spans as cpg_distribute.
// (though not perform the other checks on it)


// see: find_spans to find the entry point to the code.
// All arguments should be Vectors (i.e. arrays)
// and all should be of type REAL except r_init_n which should be an integer
// giving the expected or initial number of spans...
// as the position as in a genome we may end up with positions larger
// than 2e9.


// grow_spans; this is only for 
// I'm using xlengthgets here, but there is not actually a really good reason to use
// it as it anyway calls allocVector and then copies the data.
// This means we don't save any memory by calling xlengthgets as it does not call reAlloc
// (and cannot as the call to allocVector does not include the SEXP as an argument)

SEXP resize_spans(SEXP r_spans, int col_number){
  SEXP r_dims = PROTECT(getAttrib( r_spans, R_DimSymbol ));
  int *dims = INTEGER( r_dims );
  // we grow the column number..
  dims[1] = col_number;
  SEXP old_spans = r_spans;
  PROTECT( r_spans = xlengthgets( r_spans, dims[0] * dims[1] ) );
  UNPROTECT_PTR(old_spans);
  setAttrib( r_spans, R_DimSymbol, r_dims );
  UNPROTECT_PTR(r_dims);
  UNPROTECT_PTR(r_spans);
  return( r_spans );
}

SEXP grow_spans(SEXP r_spans, double grow_factor){
  SEXP r_dims = PROTECT(getAttrib( r_spans, R_DimSymbol ));
  int *dims = INTEGER( r_dims );
  // we grow the column number..
  dims[1] = dims[1] * grow_factor;
  SEXP old_spans = r_spans;
  PROTECT( r_spans = xlengthgets( r_spans, dims[0] * dims[1] ) );
  UNPROTECT_PTR(old_spans);
  setAttrib( r_spans, R_DimSymbol, r_dims );
  UNPROTECT_PTR(r_dims);
  return( r_spans );
}


// size_t should probably be R_xlen_t, but that doesn't make much sense if the INTEGER gives me int*
// this is a bit strange.. 
void add_span(double *scores, double *positions, size_t min_i, size_t max_i, double max_score, SEXP *r_spans, size_t *ncol, size_t *size){
  if(*size >= *ncol){
    *r_spans = grow_spans(*r_spans, 2);
    *ncol = INTEGER(getAttrib( *r_spans, R_DimSymbol ))[1];
  }
  // calculate the sum of the scores within the region
  double sum_score = 0;
  for(size_t i=min_i; i <= max_i; ++i)
    sum_score += scores[ i ];

  size_t nrow = INTEGER( getAttrib(*r_spans, R_DimSymbol) )[0];
  size_t i = (*size) * nrow;
  double *v = REAL( *r_spans ) + i;
  // we have an assumption that nrow is 6. Let's confirm that..
  if(nrow != 7){
    // note that this would probably give rise to a memory leak if we return here
    error("add_span nrow is not 6 but %d\n", nrow);
  }
  // we add 1 here as R counts from 1 not 0
  v[0] = min_i + 1;
  v[1] = max_i + 1;
  v[2] = positions[min_i];
  v[3] = positions[max_i];
  v[4] = max_score;
  v[5] = sum_score / (double)(1 + max_i - min_i);
  v[6] = (double)(1 + max_i - min_i);
  ++(*size);
}

// scores and positions 
void locate_spans(double *scores, double *positions, double *pos_breaks, size_t b_i, size_t pb_size,
		  double sep_penalty, size_t start, size_t end, SEXP *r_spans, size_t *ncol, size_t *size){
  double current_score = scores[start] > 0 ? scores[start] : 0;
  double last_score = 0;
  double max_score = current_score;
  //  double sum_score = current_score;
  ssize_t max_i = start;
  ssize_t min_i = start;
  while( b_i < pb_size && positions[ start ] > pos_breaks[ b_i ] )
    ++b_i;

  ssize_t i = start + 1;
  while(i < end){
    if( b_i < pb_size && positions[i] > pos_breaks[ b_i ] ){
      ++b_i;
      if( current_score > 0 )
	add_span( scores, positions, min_i, max_i, max_score, r_spans, ncol, size);
      current_score = 0;
      max_score = current_score;
      max_i = i;
      min_i = i;
      //      ++i;
      continue;
    }
    last_score = current_score;
    if(current_score){
      current_score = current_score + scores[i] - (sep_penalty * (positions[i] - positions[i-1]));
    }else{
      current_score = scores[i];
    }
    current_score = current_score < 0 ? 0 : current_score;
    if(current_score > max_score){
      max_score = current_score;
      max_i = i;
    }
    if(current_score > 0 && !last_score){
      min_i = i;
    }
    if(last_score && !current_score){ // last score is positive, 
      // extend the sep_penalty scores, and reset the current score, max score and so on..
      add_span(scores, positions, min_i, max_i, max_score, r_spans, ncol, size);
      // then recurse from the maximum position.. + 1; unless max_i == i
      if((i-2) > max_i)
	locate_spans(scores, positions, pos_breaks, b_i, pb_size, sep_penalty, max_i + 1, i, r_spans, ncol, size);
      // then continue on our merry way from i + 1
      ++i;
      if(i >= end)
	break;
      // current_score is set to the next score as we need to look between scores.. 
      current_score = scores[i] > 0 ? scores[i] : 0;
      max_score = current_score;
      //      sum_score = current_score;
      max_i = i;
      min_i = i;
    }
    ++i;
  }
  // Here we have run to the end; If we have a positive score we need to report that.. 
  if(current_score){
    add_span(scores, positions, min_i, max_i, max_score, r_spans, ncol, size);
    if( (end - 2) > max_i )
      locate_spans(scores, positions, pos_breaks, b_i, pb_size, sep_penalty, max_i + 1, end, r_spans, ncol, size);
  }
}

SEXP find_spans( SEXP r_scores, SEXP r_positions, SEXP r_sep_penalty, SEXP r_init_n, SEXP r_chr_offsets ){
  if(!isVector(r_scores) || !isVector(r_positions) || !isVector(r_sep_penalty) || !isVector( r_chr_offsets ))
    error("All arguments should be vectors\n");
  if(!isReal(r_scores) || !isReal(r_positions) || !isReal(r_sep_penalty) || !isReal(r_chr_offsets))
    error("All arguments should be REAL\n");

  size_t l = length(r_scores);
  if(l != length(r_positions))
    error("The scores and position vectors must be of the same size\n");

  if(length(r_sep_penalty) != 1)
    error("A single seperation penalty should be given\n");
  
  if(!isInteger(r_init_n))
    error("The last argument should be a single integer giving the expected number of spans\n");
  
  size_t pb_size = length( r_chr_offsets );
  double *chr_offsets = REAL( r_chr_offsets );

  double sep_penalty = REAL(r_sep_penalty)[0];
  double *scores = REAL(r_scores);
  double *positions = REAL(r_positions);
  
  int init_n = asInteger( r_init_n );
  if(init_n <= 0)
    error("r_init_n must be larger than 0\n");
  // we then need to create a suitable data structure for keeping the span data
  // We need to include:
  // Span starts, ends by index
  // Span starts, ends by position
  // Span scores (i.e. the max span score)
  // Span mean score (simply the average of the scores ignoring the seperation)
  
  // Start with init_n columns. Everything is REAL..
  SEXP r_spans = PROTECT(allocMatrix( REALSXP, 7, init_n));

  size_t ncol = init_n;
  size_t size = 0;
  
  locate_spans(scores, positions, chr_offsets, 0, pb_size, sep_penalty, 0, l, &r_spans, &ncol, &size);
  r_spans = PROTECT( resize_spans( r_spans, size ) );
  // here we have to set the dimensions of the vector. Should be 7 rows and size columns
  /* SEXP r_dims = getAttrib( r_spans, R_DimSymbol ); */
  /* int *dims = INTEGER( r_dims ); */
  /* dims[1] = size; */
  /* setAttrib( r_spans, R_DimSymbol, r_dims ); */
  UNPROTECT(1);
  return( r_spans );
}
