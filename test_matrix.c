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
    matrix_t *m = CreateMatrixInt(3, 4);
    ASSERT("CreateMatrixInt returns non-NULL", m != NULL);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            MatrixSetInt(m, i, j, i * 4 + j);
    ASSERT("element [2][3] == 11", MatrixGetInt(m, 2, 3) == 11);
    freedom(m);
    printf("\n");
}

static void test_create_float(void) {
    printf("-- CreateMatrixFloat --\n");
    matrix_t *m = CreateMatrixFloat(2, 2);
    ASSERT("CreateMatrixFloat returns non-NULL", m != NULL);
    MatrixSetFloat(m, 0, 0, 1.5f);
    MatrixSetFloat(m, 0, 1, 2.5f);
    MatrixSetFloat(m, 1, 0, 3.5f);
    MatrixSetFloat(m, 1, 1, 4.5f);
    ASSERT("float element [1][1] == 4.5", fabsf(MatrixGetFloat(m, 1, 1) - 4.5f) < 1e-6f);
    freedom(m);
    printf("\n");
}

static void test_transpose(void) {
    printf("-- TransMatrixInt --\n");
    matrix_t *m = CreateMatrixInt(2, 3);
    MatrixSetInt(m, 0, 0, 1); MatrixSetInt(m, 0, 1, 2); MatrixSetInt(m, 0, 2, 3);
    MatrixSetInt(m, 1, 0, 4); MatrixSetInt(m, 1, 1, 5); MatrixSetInt(m, 1, 2, 6);
    matrix_t *t = TransMatrixInt(m);
    ASSERT("transpose [0][0] == 1", MatrixGetInt(t, 0, 0) == 1);
    ASSERT("transpose [0][1] == 4", MatrixGetInt(t, 0, 1) == 4);
    ASSERT("transpose [1][0] == 2", MatrixGetInt(t, 1, 0) == 2);
    ASSERT("transpose [2][1] == 6", MatrixGetInt(t, 2, 1) == 6);
    freedom(t);
    freedom(m);
    printf("\n");
}

static void test_submatrix(void) {
    printf("-- SubMatrix --\n");
    matrix_t *m = CreateMatrixInt(3, 3);
    int v = 1;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            MatrixSetInt(m, i, j, v++);
    matrix_t *s = SubMatrix(m, 1, 1);
    ASSERT("sub[0][0] == 5", MatrixGetInt(s, 0, 0) == 5);
    ASSERT("sub[0][1] == 6", MatrixGetInt(s, 0, 1) == 6);
    ASSERT("sub[1][0] == 8", MatrixGetInt(s, 1, 0) == 8);
    ASSERT("sub[1][1] == 9", MatrixGetInt(s, 1, 1) == 9);
    freedom(s);
    freedom(m);
    printf("\n");
}

static void test_determinant(void) {
    printf("-- det --\n");
    matrix_t *m1 = CreateMatrixInt(1, 1);
    MatrixSetInt(m1, 0, 0, 7);
    ASSERT("det of [7] == 7", det(m1) == 7);
    freedom(m1);

    matrix_t *m2 = CreateMatrixInt(2, 2);
    MatrixSetInt(m2, 0, 0, 1); MatrixSetInt(m2, 0, 1, 2);
    MatrixSetInt(m2, 1, 0, 3); MatrixSetInt(m2, 1, 1, 4);
    ASSERT("det of [[1,2],[3,4]] == -2", det(m2) == -2);
    freedom(m2);

    matrix_t *m3 = CreateMatrixInt(3, 3);
    MatrixSetInt(m3, 0, 0, 1); MatrixSetInt(m3, 0, 1, 2); MatrixSetInt(m3, 0, 2, 3);
    MatrixSetInt(m3, 1, 0, 0); MatrixSetInt(m3, 1, 1, 1); MatrixSetInt(m3, 1, 2, 4);
    MatrixSetInt(m3, 2, 0, 5); MatrixSetInt(m3, 2, 1, 6); MatrixSetInt(m3, 2, 2, 0);
    ASSERT("det of 3x3 == 1", det(m3) == 1);
    freedom(m3);
    printf("\n");
}

static void test_inverse(void) {
    printf("-- Inverse --\n");
    matrix_t *m = CreateMatrixInt(2, 2);
    MatrixSetInt(m, 0, 0, 1); MatrixSetInt(m, 0, 1, 2);
    MatrixSetInt(m, 1, 0, 3); MatrixSetInt(m, 1, 1, 4);
    matrix_t *inv = Inverse(m);
    ASSERT("inv[0][0] approx -2.0", fabsf(MatrixGetFloat(inv, 0, 0) - (-2.0f)) < 0.01f);
    ASSERT("inv[0][1] approx 1.0",  fabsf(MatrixGetFloat(inv, 0, 1) - 1.0f)    < 0.01f);
    ASSERT("inv[1][0] approx 1.5",  fabsf(MatrixGetFloat(inv, 1, 0) - 1.5f)    < 0.01f);
    ASSERT("inv[1][1] approx -0.5", fabsf(MatrixGetFloat(inv, 1, 1) - (-0.5f)) < 0.01f);
    freedom(inv);
    freedom(m);
    printf("\n");
}

static void test_multiply(void) {
    printf("-- MultipliMatrix --\n");
    matrix_t *a = CreateMatrixInt(2, 3);
    MatrixSetInt(a, 0, 0, 1); MatrixSetInt(a, 0, 1, 2); MatrixSetInt(a, 0, 2, 3);
    MatrixSetInt(a, 1, 0, 4); MatrixSetInt(a, 1, 1, 5); MatrixSetInt(a, 1, 2, 6);
    matrix_t *b = CreateMatrixInt(3, 2);
    MatrixSetInt(b, 0, 0, 7);  MatrixSetInt(b, 0, 1, 8);
    MatrixSetInt(b, 1, 0, 9);  MatrixSetInt(b, 1, 1, 10);
    MatrixSetInt(b, 2, 0, 11); MatrixSetInt(b, 2, 1, 12);
    matrix_t *c = MultipliMatrix(a, b);
    ASSERT("product[0][0] == 58",  MatrixGetInt(c, 0, 0) == 58);
    ASSERT("product[0][1] == 64",  MatrixGetInt(c, 0, 1) == 64);
    ASSERT("product[1][0] == 139", MatrixGetInt(c, 1, 0) == 139);
    ASSERT("product[1][1] == 154", MatrixGetInt(c, 1, 1) == 154);
    freedom(c);
    freedom(b);
    freedom(a);
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
