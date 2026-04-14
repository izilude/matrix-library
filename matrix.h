#ifndef MATRIX_H
#define MATRIX_H

typedef struct matrix_t matrix_t;

/* Creation / Destruction */
matrix_t *CreateMatrixInt(int i, int j);
matrix_t *CreateMatrixFloat(int i, int j);
void freedom(matrix_t *m);

/* Element access */
int MatrixGetInt(matrix_t *m, int r, int c);
void MatrixSetInt(matrix_t *m, int r, int c, int val);
float MatrixGetFloat(matrix_t *m, int r, int c);
void MatrixSetFloat(matrix_t *m, int r, int c, float val);

/* I/O */
void PrintMatrixInt(matrix_t *m);
void PrintMatrixfloat(matrix_t *m);
void InputMatrixInt(matrix_t *m);

/* Operations */
matrix_t *TransMatrixInt(matrix_t *m);
matrix_t *SubMatrix(matrix_t *m, int l, int c);
int det(matrix_t *m);
matrix_t *Inverse(matrix_t *m);
matrix_t *MultipliMatrix(matrix_t *a, matrix_t *b);

#include "matrix.c"

#endif
