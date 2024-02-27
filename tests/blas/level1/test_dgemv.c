#include "test.h"

MunitResult test_dgemv_0_row(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float64_t alpha = { SB_REAL64_ONE };  // scalar alpha
    const float64_t beta = { SB_REAL64_ZERO };  // scalar beta

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float64_t* A = dvec((double[]){1.0, 2.0, 3.0, 4.0}, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float64_t* DX = dvec((double[]){1.0, 2.0}, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float64_t* DY = dvec((double[]){3.0, 4.0}, 2);

    // stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call dgemv
    dgemv(Layout, Trans, M, N, alpha, A, lda, DX, incX, beta, DY, incY);

    float64_t* RY = dvec((double[]){5.0, 11.0}, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(A);
    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_dgemv_0_col(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'C';  // Column-major order
    const char Trans = 'N';   // No transpose
    const float64_t alpha = { SB_REAL64_ONE };  // Dcalar alpha
    const float64_t beta = { SB_REAL64_ZERO };  // Dcalar beta

    // 2x2 matrix A (column-major order)
    //   [1.0 3.0]
    //   [2.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float64_t* A = dvec((double[]){1.0, 3.0, 2.0, 4.0}, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float64_t* DX = dvec((double[]){1.0, 2.0}, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float64_t* DY = dvec((double[]){3.0, 4.0}, 2);

    // Dtride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = M;

    // Call dgemv
    dgemv(Layout, Trans, M, N, alpha, A, lda, DX, incX, beta, DY, incY);

    float64_t* RY = dvec((double[]){5.0, 11.0}, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(A);
    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_dgemv_12345(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float64_t alpha = {*(uint64_t*)&(double){2.5f}};
    const float64_t beta = {*(uint64_t*)&(double){-1.0}};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 4;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float64_t* A = dvec((double[]){ 2.0, -1.0,  3.0,  0.0,  4.0,
                                   1.0,  0.0,  2.0, -1.0,  3.0,
                                   4.0,  5.0,  1.0,  2.0,  0.0,
                                   3.0,  2.0,  4.0,  2.0,  5.0},
                        M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float64_t* DY = dvec((double[]){-3.0, 4.0, -5.0, 6.0f}, 4);

    // Dtride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call dgemv
    dgemv(Layout, Trans, M, N, alpha, A, lda, DX, incX, beta, DY, incY);

    float64_t* RY = dvec((double[]){75.5, 41.0, 67.5, 124.0}, 4);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(A);
    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_dgemv_stride(const MunitParameter params[],
                              void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float64_t alpha = {*(uint64_t*)&(double){2.5f}};
    const float64_t beta = {*(uint64_t*)&(double){-1.0}};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 7;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float64_t* A = dvec((double[]){ 2.0, -1.0,  3.0,  0.0,  4.0,
                                   0.0,  0.0,  0.0,  0.0,  0.0,
                                   1.0,  0.0,  2.0, -1.0,  3.0,
                                   0.0,  0.0,  0.0,  0.0,  0.0,
                                   4.0,  5.0,  1.0,  2.0,  0.0,
                                   0.0,  0.0,  0.0,  0.0,  0.0,
                                   3.0,  2.0,  4.0,  2.0,  5.0},
                       M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float64_t* DY = dvec((double[]){-3.0, 0.0, 4.0, 0.0, -5.0, 0.0, 6.0f}, 7);

    // Dtride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 2;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call dgemv
    dgemv(Layout, Trans, 4, 5, alpha, A, lda, DX, incX, beta, DY, incY);

    float64_t* RY = dvec((double[]){75.5, 0.0, 41.0, 0.0, 67.5, 0.0, 124.0}, 7);

    for (uint64_t i = 0; i < 7; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(A);
    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}
