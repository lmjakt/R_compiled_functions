#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>

// provide a function to read in a set of fasta and find instances of specific
// patterns...

// read_one_sequence based on code from Gos Micklem and Tim Cutts cpg_distribute.c
// and hence presumably falls under the GPL version 2.1 or later licensing.

void rm_eol(const char *word){
  char *p;
  for(p=(char*)word; *p != '\n'; p++);
  *p = 0;
}

char* read_one_sequence(FILE *f, char* title, size_t title_size, char *seq_buffer, size_t *seq_buffer_size, size_t *len )
{
  char *p;
  char buf[1024];
  int n;
  size_t last_offset = 0;

  *seq_buffer = '\0'; *len = 0;
  *title = '\0';

  while (!feof(f)) {
    last_offset = ftell( f );
    p = fgets(buf, sizeof(buf), f);
    if (p) {
      if (*p == '>') {
	/* We have found a comment line, so we can process */
	if (*seq_buffer != '\0') {
	  seq_buffer[*len] = '\0';
	  fseek( f, last_offset, SEEK_SET );
	  return( seq_buffer );
	}
	strncpy(title, ++p, title_size);
	rm_eol( title );
	continue;
      }
      n = strlen(p)-1;
      if ((*len+n+3) > *seq_buffer_size) {
	*seq_buffer_size <<= 1;
	seq_buffer = (char *)realloc(seq_buffer, *seq_buffer_size);
      }
      if (!seq_buffer) {
	fprintf(stderr, "Out of memory\n");
	return(seq_buffer);
      }
      memcpy(&seq_buffer[*len], buf, n);
      *len += n;
    }
  }
  
  if (*seq_buffer != '\0') {
    seq_buffer[*len] = '\0';
    return( seq_buffer );
  }
  return(seq_buffer);
}

SEXP resize_vector(SEXP v, size_t new_size){
  SEXP old_v = v;
  PROTECT( v = xlengthgets( v, new_size ) );
  UNPROTECT_PTR( old_v );
  UNPROTECT_PTR( v ); // we unprotect in case the return value is protected by being added to a vector.
  return( v );
}

SEXP scan_sequence(const char *seq, const char *word, size_t init_length){
  SEXP positions = PROTECT(allocVector( INTSXP, init_length ) );
  
  size_t pos_l = length( positions );
  int *pos = INTEGER( positions );
  size_t l = strlen( word );
  size_t j = 0;
  for(size_t i = 0; seq[i] != 0; ++i){
    if(!strncmp( word, seq + i, l )){
      if(j >= pos_l){
	positions = PROTECT(resize_vector( positions, pos_l * 2 ));
	pos_l = length( positions );
	pos = INTEGER( positions );
      }
      pos[ j++ ] = i;
    }
  }
  return( resize_vector( positions, j ) );
}


SEXP find_words(SEXP r_fname, SEXP r_word){
  // both fname an word should be a single character string
  if( !isString(r_fname) || !isString(r_word) )
    error("both arguments should be character vectors\n");
  if(length(r_fname) <= 0 || length(r_word) <= 0)
    error("both arguments must contain at least a single word\n");

  const char *fname = (const char*) CHAR( STRING_ELT( r_fname, 0 ) );
  const char *word = (const char*) CHAR( STRING_ELT( r_word, 0 ) );

  FILE *file;
  if(!(file = fopen( fname, "r" )))
    error("Unable to open file");
  
  // allocate some data structures that we will use.. 
  size_t seq_buffer_size = 8192;
  char *seq_buffer = (char*)malloc( seq_buffer_size * sizeof(char) );
  size_t title_size = 1024;
  char *title = (char*)malloc( title_size * sizeof(char) );
  size_t len = 0;

  size_t vec_length=50;
  size_t match_init_length = 1000;
  SEXP seq_names = PROTECT( allocVector( STRSXP, vec_length ) );
  SEXP seq_lengths = PROTECT( allocVector( INTSXP, vec_length ) );
  SEXP seq_matches = PROTECT( allocVector( VECSXP, vec_length ) );

  int *seq_l = INTEGER( seq_lengths );
  size_t seq_no = 0;

  while(1){
    seq_buffer = read_one_sequence( file, title, title_size, seq_buffer, &seq_buffer_size, &len );
    if(len == 0)
      break;
    
    if(seq_no >= vec_length){
      vec_length = 2 * vec_length;
      seq_names = PROTECT( resize_vector( seq_names, vec_length ) );
      seq_lengths = PROTECT( resize_vector( seq_lengths, vec_length ) );
      seq_matches = PROTECT( resize_vector( seq_matches, vec_length ) );
      seq_l = INTEGER( seq_lengths );
    }
    SET_STRING_ELT( seq_names, seq_no, mkChar( title ) );
    seq_l[ seq_no ] = len;
    SET_VECTOR_ELT( seq_matches, seq_no, scan_sequence( seq_buffer, word, match_init_length) );

    seq_no++;
  }
  
  // unfortunately this uses malloc() and copy, rather than realloc
  // so ends up being expensive.. 
  seq_names = PROTECT( resize_vector( seq_names, seq_no ) );
  seq_lengths = PROTECT( resize_vector( seq_lengths, seq_no ) );
  seq_matches = PROTECT( resize_vector( seq_matches, seq_no ) );

  free( seq_buffer );
  free( title );

  SEXP ret_data = PROTECT( allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( ret_data, 0, seq_names );
  SET_VECTOR_ELT( ret_data, 1, seq_lengths );
  SET_VECTOR_ELT( ret_data, 2, seq_matches );
  
  UNPROTECT(4);
  return( ret_data );

}

