#ifndef MATRIX_C
#define MATRIX_C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct matrix_t {
    void *data;
    int rows;
    int cols;
    int elem_type;        /* 0 = int, 1 = float */
};

/* ================================================================
 * Creation / Destruction
 * ================================================================ */

matrix_t *CreateMatrixInt(int i, int j)
{
    int **d = (int **)malloc(i * sizeof(int *));
    for (int r = 0; r < i; r++)
        d[r] = (int *)calloc(j, sizeof(int));

    matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
    m->data = d;
    m->rows = i;
    m->cols = j;
    m->elem_type = 0;
    return m;
}

matrix_t *CreateMatrixFloat(int i, int j)
{
    float **d = (float **)malloc(i * sizeof(float *));
    for (int r = 0; r < i; r++)
        d[r] = (float *)calloc(j, sizeof(float));

    matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
    m->data = d;
    m->rows = i;
    m->cols = j;
    m->elem_type = 1;
    return m;
}

void freedom(matrix_t *m)
{
    if (!m) return;
    if (m->elem_type == 0)
    {
        int **d = (int **)m->data;
        for (int i = 0; i < m->rows; i++)
            free(d[i]);
        free(d);
    }
    else
    {
        float **d = (float **)m->data;
        for (int i = 0; i < m->rows; i++)
            free(d[i]);
        free(d);
    }
    free(m);
}

/* ================================================================
 * Element access
 * ================================================================ */

int MatrixGetInt(matrix_t *m, int r, int c)
{
    return ((int **)m->data)[r][c];
}

void MatrixSetInt(matrix_t *m, int r, int c, int val)
{
    ((int **)m->data)[r][c] = val;
}

float MatrixGetFloat(matrix_t *m, int r, int c)
{
    return ((float **)m->data)[r][c];
}

void MatrixSetFloat(matrix_t *m, int r, int c, float val)
{
    ((float **)m->data)[r][c] = val;
}

/* ================================================================
 * I/O
 * ================================================================ */

void InputMatrixInt(matrix_t *m)
{
    int **d = (int **)m->data;
    for (int i = 0; i < m->rows; i++)
    {
        printf("//Digite os %d elementos da linha %d\n", m->cols, (i + 1));
        for (int j = 0; j < m->cols; j++)
        {
            scanf("%d", &d[i][j]);
        }
    }
}

void PrintMatrixInt(matrix_t *m)
{
    int **d = (int **)m->data;
    for (int i = 0; i < m->rows; i++)
    {
        for (int j = 0; j < m->cols; j++)
        {
            printf("%5d", d[i][j]);
        }
        printf("\n");
    }
}

void PrintMatrixfloat(matrix_t *m)
{
    float **d = (float **)m->data;
    for (int i = 0; i < m->rows; i++)
    {
        for (int j = 0; j < m->cols; j++)
        {
            printf("%10.3g", d[i][j]);
        }
        printf("\n");
    }
}

/* ================================================================
 * Operations
 * ================================================================ */

matrix_t *TransMatrixInt(matrix_t *m)
{
    int **d = (int **)m->data;
    matrix_t *t = CreateMatrixInt(m->cols, m->rows);
    int **td = (int **)t->data;
    for (int i = 0; i < m->cols; i++)
    {
        for (int j = 0; j < m->rows; j++)
        {
            td[i][j] = d[j][i];
        }
    }
    return t;
}

matrix_t *SubMatrix(matrix_t *m, int l, int c)
{
    int o = m->rows;
    int **d = (int **)m->data;
    matrix_t *sub = CreateMatrixInt(o - 1, o - 1);
    int **sd = (int **)sub->data;
    for (int i = 0, w = 0; i < o; i++)
    {
        if (i == (l - 1))
            continue;
        else
        {
            w++;
            for (int j = 0, z = 0; j < o; j++)
            {
                if (j == (c - 1))
                    continue;
                else
                {
                    z++;
                    sd[w - 1][z - 1] = d[i][j];
                }
            }
        }
    }
    return sub;
}

int det(matrix_t *m)
{
    int o = m->rows;
    int **d = (int **)m->data;
    int r = 0;
    if (o == 1)
    {
        r = d[0][0];
        return r;
    }
    else if (o == 2)
    {
        r = ((d[0][0] * d[1][1]) - (d[0][1] * d[1][0]));
        return r;
    }
    else
    {
        for (int i = 0, j = 0; j < o; j++)
        {
            r += (d[i][j] * (pow((-1), (i + j + 2))) * det(SubMatrix(m, 1, j + 1)));
        }
        return r;
    }
}

matrix_t *Inverse(matrix_t *m)
{
    int o = m->rows;
    float determinante = (float)det(m);
    matrix_t *MatrizCof = CreateMatrixInt(o, o);
    for (int i = 0; i < o; i++)
    {
        for (int j = 0; j < o; j++)
        {
            MatrixSetInt(MatrizCof, i, j,
                (int)((pow((-1), (i + j + 2))) * det(SubMatrix(m, i + 1, j + 1))));
        }
    }
    matrix_t *adjunta = TransMatrixInt(MatrizCof);
    matrix_t *inversa = CreateMatrixFloat(o, o);
    for (int i = 0; i < o; i++)
    {
        for (int j = 0; j < o; j++)
        {
            MatrixSetFloat(inversa, i, j,
                (float)(1.00 / determinante) * (float)MatrixGetInt(adjunta, i, j));
        }
    }
    freedom(adjunta);
    freedom(MatrizCof);
    return inversa;
}

matrix_t *MultipliMatrix(matrix_t *a, matrix_t *b)
{
    int la = a->rows, ca = a->cols, cb = b->cols;
    int **da = (int **)a->data;
    int **db = (int **)b->data;
    matrix_t *result = CreateMatrixInt(la, cb);
    int **dr = (int **)result->data;
    for (int i = 0; i < la; i++)
    {
        for (int j = 0; j < cb; j++)
        {
            dr[i][j] = 0;
            for (int w = 0; w < ca; w++)
            {
                dr[i][j] += (da[i][w] * db[w][j]);
            }
        }
    }
    return result;
}

#endif
