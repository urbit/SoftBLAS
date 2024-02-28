#include "test.h"
#include <stdio.h>

MunitResult test_dgemm_0_row(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float64_t alpha = { SB_REAL64_ONE };  // Dcalar alpha
    const float64_t beta = { SB_REAL64_ZERO };  // Dcalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float64_t* A = dvec((double[]){1.0, 2.0, 3.0, 4.0}, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float64_t* B = dvec((double[]){5.0, 6.0, 7.0, 8.0}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float64_t* C = dvec((double[]){0.0, 0.0, 0.0, 0.0}, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call dgemm
    dgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float64_t* R = dvec((double[]){19.0, 22.0, 43.0, 50.0}, M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_dgemm_0_col(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'T';  // Transpose
    const char transB = 'T';  // Transpose
    const float64_t alpha = { SB_REAL64_ONE };  // Dcalar alpha
    const float64_t beta = { SB_REAL64_ZERO };  // Dcalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float64_t* A = dvec((double[]){1.0, 3.0, 2.0, 4.0}, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float64_t* B = dvec((double[]){5.0, 7.0, 6.0, 8.0}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float64_t* C = dvec((double[]){0.0, 0.0, 0.0, 0.0}, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call dgemm
    dgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float64_t* R = dvec((double[]){19.0, 22.0, 43.0, 50.0}, M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_dgemm_3x2x1_0(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';  // No transpose
    const char transB = 'N';  // No transpose
    const float64_t alpha = {*(uint64_t*)&(double){1.0}};
    const float64_t beta = {*(uint64_t*)&(double){0.0}};

    const uint64_t M = 3;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 1;     // Number of columns in B

    // 3x2 matrix A (row-major order)
    float64_t* A = dvec((double[]){ 1.0,  2.0,
                                   3.0,  4.0,
                                   5.0,  6.0},
                        M*N);

    // 2x1 matrix B (row-major order)
    float64_t* B = dvec((double[]){ 1.0,
                                   2.0},
                        N*P);

    // 3x1 matrix C (row-major order)
    float64_t* C = dvec((double[]){ 0.0,
                                   0.0,
                                   0.0},
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call dgemm
    dgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float64_t* R = dvec((double[]){ 5.0,
                                   11.0,
                                   17.0},
                        M*P);


    for (uint64_t i = 0; i < 3; i++) {
        assert_ullong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_dgemm_5x4x3(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';  // No transpose
    const char transB = 'N';  // No transpose
    const float64_t alpha = {*(uint64_t*)&(double){2.5}};
    const float64_t beta = {*(uint64_t*)&(double){-1.0}};

    const uint64_t M = 5;     // Number of rows in A
    const uint64_t N = 4;     // Number of columns in A
    const uint64_t P = 3;     // Number of columns in B

    // 5x4 matrix A (row-major order)
    float64_t* A = dvec((double[]){ 2.0, -1.0,  3.0,  0.0,
                                   4.0,  1.0,  0.0,  2.0,
                                  -1.0,  3.0,  4.0,  5.0,
                                   1.0,  2.0,  0.0,  3.0,
                                   2.0,  4.0,  2.0,  5.0},
                        M*N);

    // 4x3 matrix B (row-major order)
    float64_t* B = dvec((double[]){ 6.0, -2.0,  8.0,
                                   4.0,  1.0,  7.0,
                                  -1.0,  2.0,  4.0,
                                   5.0,  4.0,  2.0},
                        N*P);

    // 2x2 matrix C (row-major order)
    float64_t* C = dvec((double[]){ 1.0,  2.0,  3.0,
                                   4.0,  5.0,  6.0,
                                  -1.0, -2.0, -3.0,
                                  -4.0, -5.0, -6.0,
                                  -2.0, -3.0, -4.0},
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call dgemm
    dgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float64_t* R = dvec((double[]){11.5,  0.5,  49.5,
                                  91.0, -2.5, 101.5,
                                  68.5, 84.5, 100.5,
                                  76.5, 35.0,  76.0,
                                 129.5, 63.0, 159.0},
                        M*P);

    for (uint64_t i = 0; i < M*P; i++) {
        assert_ullong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}
