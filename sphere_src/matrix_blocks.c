#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "matrix_blocks.h"

/*
 * Matrices are read into a list of chunks of the following number of
 * floating point elements.
 */
#define MTX_BLK_SIZE 1024

/*
 * Objects of the following type are used to record a single chunk
 * of a matrix of floating point numbers.
 */
typedef struct MatrixBlk MatrixBlk;
struct MatrixBlk {
  MatrixBlk *next;           /* A pointer to the next chunk of the matrix */
  float array[MTX_BLK_SIZE]; /* An array containing one part of a matrix */
  int n;                     /* The number of elements currently in array[] */
};

/*
 * Objects of the following type are used to record lists of blocks of
 * matrices, as their constituent numbers are are read in from a file
 * of unknown size.
 */
struct MatrixBlocks {
  MatrixBlk *head;        /* The head of the list of sequential chunks */
  MatrixBlk *tail;        /* The tail of the list of sequential chunks */
  int nrow;               /* The number of rows in the matrix */
  int ncol;               /* The number of columns in the matrix */
};

static int mb_read_file(MatrixBlocks *mb, const char *file, int nskip);
static int mb_cant_read_file(FILE *fp);

/*.......................................................................
 * Create a new MatrixBlocks object, filled with a matrix read from a
 * file containing an unknown number of rows, of an unknown number of
 * columns.
 *
 * Input:
 *  file      const char *   The name of the file to read the matrix
 *                           from.
 *  nskip            int     The number of header lines to skip before
 *                           attempting to read columns from 'file'.
 * Output:
 *  return  MatrixBlocks *  The new object, or NULL on error.
 */
MatrixBlocks *new_MatrixBlocks(const char *file, int nskip)
{
  MatrixBlocks *mb;  /* The object to be returned */
/*
 * Check the arguments.
 */
  if(!file) {
    lprintf(stderr, "new_MatrixBlocks: NULL filename string.\n");
    return NULL;
  };
/*
 * Allocate the container.
 */
  mb = malloc(sizeof(MatrixBlocks));
  if(!mb) {
    lprintf(stderr, "new_MatrixBlocks: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_MatrixBlocks().
 */
  mb->head = NULL;
  mb->tail = NULL;
  mb->nrow = 0;
  mb->ncol = 0;
/*
 * Allocate the first node of the list of blocks.
 */
  mb->head = mb->tail = malloc(sizeof(MatrixBlk));
  if(!mb->head) {
    lprintf(stderr, "new_MatrixBlocks: Insufficient memory.\n");
    return del_MatrixBlocks(mb);
  };
/*
 * Mark the first block as currently being empty.
 */
  mb->head->n = 0;
  mb->head->next = NULL;
/*
 * Read the file.
 */
  if(mb_read_file(mb, file, nskip))
    return del_MatrixBlocks(mb);
  return mb;
}

/*.......................................................................
 * Delete a MatrixBlocks object.
 *
 * Input:
 *  mb     MatrixBlocks *  The object to be deleted.
 * Output:
 *  return MatrixBlocks *  The deleted object (always NULL).
 */
MatrixBlocks *del_MatrixBlocks(MatrixBlocks *mb)
{
  if(mb) {
/*
 * Free the nodes of the list of blocks.
 */
    while(mb->head) {
      MatrixBlk *node = mb->head;
      mb->head = node->next;
      free(node);
    };
    mb->head = NULL;
    mb->tail = NULL;
/*
 * Discard the container.
 */
    free(mb);
  };
  return NULL;
}

/*.......................................................................
 * This is a private function of new_MatrixBlocks(), used to read a file
 * that contains an unknown number of lines and columns. Note that the
 * number of columns must not change from one line to the next.
 *
 * Input:
 *  mb    MatrixBlocks *   The container into which to read the numbers.
 *  nskip          int     The number of lines to skip at the start of
 *                         the file.
 *  file    const char *   The file to read from.
 * Output:
 *  return         int     0 - OK.
 *                         1 - Error.
 */
static int mb_read_file(MatrixBlocks *mb, const char *file, int nskip)
{
  FILE *fp;   /* The stream of the open file */
  int c='\0'; /* A character read from the file */
  int ncol;   /* The number of columns read from the current line */
/*
 * Attempt to open the file.
 */
  fp = fopen(file, "r");
  if(!fp) {
    lprintf(stderr, "Unable to open file: %s\n", file);
    return 1;
  };
/*
 * Skip nskip lines.
 */
  while(nskip > 0 && c!=EOF) {
    do {
      c = fgetc(fp);
    } while(c != '\n' && c != EOF);
    nskip--;
  };
/*
 * Read numbers, while checking that the number of numbers per line is
 * consistent. Record the numbers in the list of matrix blocks, allocating
 * more blocks, as needed.
 */
  ncol = 0;
  while(!feof(fp)) {
    double d;         /* The latest number read from the file */
/*
 * Skip leading spaces and tabs.
 */
    do { c = fgetc(fp); } while(c==' ' || c=='\t');
/*
 * If we hit a newline character, this represents the end of a new line.
 * If this is the first line, the number of columns read from this line
 * sets the number expected from all subsequent lines. Check that this
 * is true on subsequent lines. Ignore completely empty lines.
 */
    while((c=='\n' || c==EOF) && ncol!=0) {
      if(mb->nrow <= 0) {
	mb->ncol = ncol;
      } else if(mb->ncol != ncol) {
	lprintf(stderr, "Inconsistent number of columns in file: %s\n", file);
	return mb_cant_read_file(fp);
      };
/*
 * Start a new row.
 */
      mb->nrow++;
      ncol = 0;
/*
 * Skip leading spaces and tabs on the new line.
 */
      if(c!=EOF) do { c = fgetc(fp); } while(c==' ' || c=='\t');
/*
 * Did we hit the end of the file?
 */
      if(c == EOF)
	return 0;
    };
/*
 * We will have read the first character of the next number, so put
 * it back for scanf() to read.
 */
    ungetc(c, fp);
/*
 * Attempt to read the next number.
 */
    if(fscanf(fp, "%lf", &d) != 1) {
      lprintf(stderr, "Missing number in column %d, row %d of matrix in: %s\n",
	      ncol+1, mb->nrow+1, file);
      return mb_cant_read_file(fp);
    };
/*
 * If there isn't sufficient space to record the number in the latest
 * block, allocate a new one, and append it to the list.
 */
    if(mb->tail->n >= MTX_BLK_SIZE) {
      MatrixBlk *node = (MatrixBlk *) malloc(sizeof(MatrixBlk));
      if(!node) {
	lprintf(stderr, "new_MatrixBlocks: Insufficient memory.\n");
	return 1;
      };
      node->next = NULL;
      node->n = 0;
      mb->tail->next = node;
      mb->tail = node;
    };
/*
 * Record the latest number.
 */
    mb->tail->array[mb->tail->n++] = d;
/*
 * We just read another column of the current row.
 */
    ncol++;
  };
/*
 * Was the end of the file reached before any numbers were read?
 */
  if(mb->nrow==0) {
    lprintf(stderr, "No numbers were read before the end of file: %s\n", file);
    return mb_cant_read_file(fp);
  };
/*
 * Close the file.
 */
  fclose(fp);
  return 0;
}

/*.......................................................................
 * This is a private error return function of mb_read_file(), used to
 * cleanup after an error, and return the error code of its caller.
 */
static int mb_cant_read_file(FILE *fp)
{
  if(fp)
    (void) fclose(fp);
  return 1;
}

/*.......................................................................
 * Copy the contents of a list of matrix blocks into a given float array.
 *
 * Input:
 *  mb     MatrixBlocks *   The container of the list of chunks of the
 *                          matrix.
 *  matrx         float *   The array into which to copy the matrix.
 *  dim             int     The allocated dimension of matrx[]. If the
 *                          number of elements in the matrix exceeds this
 *                          number, the copy will be silently truncated.
 * Output:
 *  return          int     0 - OK.
 *                          1 - Error.
 */
int mb_to_float_array(MatrixBlocks *mb, float *matrx, int dim)
{
  int offset;     /* The index of matrx[] into which to copy the next element */
  MatrixBlk *blk; /* One block of the list of matrix blocks */
  int i;
/*
 * Copy the matrix values into the return matrix.
 */
  for(offset=0,blk=mb->head; blk; blk=blk->next) {
    int n = offset + blk->n <= dim ? blk->n : (dim - offset);
    for(i=0; i<n; i++)
      matrx[offset++] = blk->array[i];
  };
  return 0;
}

/*.......................................................................
 * Return the dimensions of the matrix contained in a given list of matrix
 * blocks read from a file.
 *
 * Input:
 *  mb     MatrixBlocks *   The container of the matrix.
 * Input/Output:
 *  nrow            int *   On output the number of rows in the matrix
 *                          will be assigned to *nrow.
 *  ncol            int *   On output the number of columns in the matrix
 *                          will be assigned to *ncol.
 */
void mb_matrix_size(MatrixBlocks *mb, int *nrow, int *ncol)
{
  if(nrow)
    *nrow = mb ? mb->nrow : 0;
  if(ncol)
    *ncol = mb ? mb->ncol : 0;
}
