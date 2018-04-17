#ifndef matrix_blocks_h
#define matrix_blocks_h

typedef struct MatrixBlocks MatrixBlocks;

MatrixBlocks *new_MatrixBlocks(const char *file, int nskip);
MatrixBlocks *del_MatrixBlocks(MatrixBlocks *mb);

void mb_matrix_size(MatrixBlocks *mb, int *nrow, int *ncol);

int mb_to_float_array(MatrixBlocks *mb, float *matrx, int dim);

#endif
