#include "test.h"
#include <stdio.h>

MunitResult test_qgemv_0_row(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float128_t alpha = { SB_REAL128L_ONE, SB_REAL128U_ONE };
    const float128_t beta = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000}},
        2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        2);

    // stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call qgemv
    qgemv(Layout, Trans, M, N, alpha, A, lda, QX, incX, beta, QY, incY);

    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002600000000000, .lo = 0x0000000000000000}},
        2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(A);
    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qgemv_0_col(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'C';  // Column-major order
    const char Trans = 'N';   // No transpose
    const float128_t alpha = { SB_REAL128L_ONE, SB_REAL128U_ONE };
    const float128_t beta = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000}},
        2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        2);

    // stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call qgemv
    qgemv(Layout, Trans, M, N, alpha, A, lda, QX, incX, beta, QY, incY);

    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002600000000000, .lo = 0x0000000000000000}},
        2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(A);
    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qgemv_12345(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float128_t alpha = { 0x0000000000000000, 0x4000400000000000 };
    const float128_t beta = { 0x0000000000000000, 0xbfff000000000000 };

    // 4x5 matrix A (row-major order)
    const uint64_t M = 4;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        M*N);

    // Vector X (column vector)
    //   [1.0, 2.0, 3.0, 4.0, 5.0]'
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0xc000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000}},
        4);

    // stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call qgemv
    qgemv(Layout, Trans, M, N, alpha, A, lda, QX, incX, beta, QY, incY);

    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x40052e0000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004480000000000, .lo = 0x0000000000000000},
            {.hi = 0x40050e0000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005f00000000000, .lo = 0x0000000000000000}},
        4);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(A);
    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qgemv_stride(const MunitParameter params[],
                              void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float128_t alpha = { 0x0000000000000000, 0x4000400000000000 };
    const float128_t beta = { 0x0000000000000000, 0xbfff000000000000 };

    // 4x5 matrix A (row-major order)
    const uint64_t M = 7;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},

            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},

            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},

            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},

            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},

            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},

            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        M*N);

    // Vector X (column vector)
    //   [1.0, 2.0, 3.0, 4.0, 5.0]'
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    // Vector Y (column vector)
    //   [-3.0 0.0 4.0 0.0 -5.0 0.0 6.0]'
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0xc000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000}},
        7);

    // stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 2;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call qgemv
    qgemv(Layout, Trans, M, N, alpha, A, lda, QX, incX, beta, QY, incY);

    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x40052e0000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004480000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x40050e0000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005f00000000000, .lo = 0x0000000000000000}},
        7);

    for (uint64_t i = 0; i < 7; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(A);
    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}
