#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

static int tests_run = 0;
static int tests_passed = 0;

#define ASSERT(msg, expr) do { \
    tests_run++; \
    if (expr) { tests_passed++; printf("  PASS: %s\n", msg); } \
    else { printf("  FAIL: %s\n", msg); } \
} while (0)

static void test_create_and_free(void) {
    printf("-- CreateMatrixInt / freedom --\n");
    int **m = CreateMatrixInt(3, 4);
    ASSERT("CreateMatrixInt returns non-NULL", m != NULL);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            m[i][j] = i * 4 + j;
    ASSERT("element [2][3] == 11", m[2][3] == 11);
    freedom(m, 3);
    printf("\n");
}

static void test_create_float(void) {
    printf("-- CreateMatrixFloat --\n");
    float **m = CreateMatrixFloat(2, 2);
    ASSERT("CreateMatrixFloat returns non-NULL", m != NULL);
    m[0][0] = 1.5f; m[0][1] = 2.5f;
    m[1][0] = 3.5f; m[1][1] = 4.5f;
    ASSERT("float element [1][1] == 4.5", fabsf(m[1][1] - 4.5f) < 1e-6f);
    for (int i = 0; i < 2; i++) free(m[i]);
    free(m);
    printf("\n");
}

static void test_transpose(void) {
    printf("-- TransMatrixInt --\n");
    int **m = CreateMatrixInt(2, 3);
    /* [ 1 2 3 ]
       [ 4 5 6 ] */
    m[0][0]=1; m[0][1]=2; m[0][2]=3;
    m[1][0]=4; m[1][1]=5; m[1][2]=6;
    int **t = TransMatrixInt(m, 3, 2);
    /* Expected (3x2):
       [ 1 4 ]
       [ 2 5 ]
       [ 3 6 ] */
    ASSERT("transpose [0][0] == 1", t[0][0] == 1);
    ASSERT("transpose [0][1] == 4", t[0][1] == 4);
    ASSERT("transpose [1][0] == 2", t[1][0] == 2);
    ASSERT("transpose [2][1] == 6", t[2][1] == 6);
    freedom(t, 3);
    freedom(m, 2);
    printf("\n");
}

static void test_submatrix(void) {
    printf("-- SubMatrix --\n");
    int **m = CreateMatrixInt(3, 3);
    /* [ 1 2 3 ]
       [ 4 5 6 ]
       [ 7 8 9 ] */
    int v = 1;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m[i][j] = v++;
    /* Remove row 1, col 1 (1-indexed) -> exclude row 0, col 0 */
    int **s = SubMatrix(m, 3, 1, 1);
    /* Expected (2x2):
       [ 5 6 ]
       [ 8 9 ] */
    ASSERT("sub[0][0] == 5", s[0][0] == 5);
    ASSERT("sub[0][1] == 6", s[0][1] == 6);
    ASSERT("sub[1][0] == 8", s[1][0] == 8);
    ASSERT("sub[1][1] == 9", s[1][1] == 9);
    freedom(s, 2);
    freedom(m, 3);
    printf("\n");
}

static void test_determinant(void) {
    printf("-- det --\n");
    /* 1x1 */
    int **m1 = CreateMatrixInt(1, 1);
    m1[0][0] = 7;
    ASSERT("det of [7] == 7", det(m1, 1) == 7);
    freedom(m1, 1);

    /* 2x2: det([[1,2],[3,4]]) = 1*4 - 2*3 = -2 */
    int **m2 = CreateMatrixInt(2, 2);
    m2[0][0]=1; m2[0][1]=2;
    m2[1][0]=3; m2[1][1]=4;
    ASSERT("det of [[1,2],[3,4]] == -2", det(m2, 2) == -2);
    freedom(m2, 2);

    /* 3x3: det([[1,2,3],[0,1,4],[5,6,0]]) = 1 */
    int **m3 = CreateMatrixInt(3, 3);
    m3[0][0]=1; m3[0][1]=2; m3[0][2]=3;
    m3[1][0]=0; m3[1][1]=1; m3[1][2]=4;
    m3[2][0]=5; m3[2][1]=6; m3[2][2]=0;
    ASSERT("det of 3x3 == 1", det(m3, 3) == 1);
    freedom(m3, 3);
    printf("\n");
}

static void test_inverse(void) {
    printf("-- Inverse --\n");
    /* Inverse of [[1,2],[3,4]] with det=-2:
       adj = [[4,-2],[-3,1]]  (cofactor transposed)
       inv = adj / det = [[-2, 1],[1.5, -0.5]] */
    int **m = CreateMatrixInt(2, 2);
    m[0][0]=1; m[0][1]=2;
    m[1][0]=3; m[1][1]=4;
    float **inv = Inverse(m, 2);
    ASSERT("inv[0][0] approx -2.0", fabsf(inv[0][0] - (-2.0f)) < 0.01f);
    ASSERT("inv[0][1] approx 1.0",  fabsf(inv[0][1] - 1.0f)    < 0.01f);
    ASSERT("inv[1][0] approx 1.5",  fabsf(inv[1][0] - 1.5f)    < 0.01f);
    ASSERT("inv[1][1] approx -0.5", fabsf(inv[1][1] - (-0.5f)) < 0.01f);
    for (int i = 0; i < 2; i++) free(inv[i]);
    free(inv);
    freedom(m, 2);
    printf("\n");
}

static void test_multiply(void) {
    printf("-- MultipliMatrix --\n");
    /* A (2x3) * B (3x2) = C (2x2) */
    int **a = CreateMatrixInt(2, 3);
    a[0][0]=1; a[0][1]=2; a[0][2]=3;
    a[1][0]=4; a[1][1]=5; a[1][2]=6;
    int **b = CreateMatrixInt(3, 2);
    b[0][0]=7;  b[0][1]=8;
    b[1][0]=9;  b[1][1]=10;
    b[2][0]=11; b[2][1]=12;
    int **c = MultipliMatrix(a, 2, 3, b, 3, 2);
    /* C = [[58,64],[139,154]] */
    ASSERT("product[0][0] == 58",  c[0][0] == 58);
    ASSERT("product[0][1] == 64",  c[0][1] == 64);
    ASSERT("product[1][0] == 139", c[1][0] == 139);
    ASSERT("product[1][1] == 154", c[1][1] == 154);
    freedom(c, 2);
    freedom(b, 3);
    freedom(a, 2);
    printf("\n");
}

int main(void) {
    printf("=== Matrix Library Tests ===\n\n");
    test_create_and_free();
    test_create_float();
    test_transpose();
    test_submatrix();
    test_determinant();
    test_inverse();
    test_multiply();
    printf("=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
