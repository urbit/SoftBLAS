#include "test.h"
#include <stdio.h>

MunitResult test_qgemm_0_row(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float128_t alpha = { SB_REAL128L_ONE, SB_REAL128U_ONE };
    const float128_t beta = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float128_t* B = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001c00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002000000000000, .lo = 0x0000000000000000}},
        N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float128_t* C = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call qgemm
    qgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float128_t* R = qvec((float128_pair_t[]){
            {.hi = 0x4003300000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003600000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004580000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004900000000000, .lo = 0x0000000000000000}},
        M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(C[i].v[0], ==, R[i].v[0]);
        assert_ullong(C[i].v[1], ==, R[i].v[1]);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_qgemm_0_col(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'T';   // No transpose
    const char transB = 'T';   // No transpose
    const float128_t alpha = { SB_REAL128L_ONE, SB_REAL128U_ONE };
    const float128_t beta = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000}},
        M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float128_t* B = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001c00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002000000000000, .lo = 0x0000000000000000}},
        N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float128_t* C = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call qgemm
    qgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float128_t* R = qvec((float128_pair_t[]){
            {.hi = 0x4003300000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003600000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004580000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004900000000000, .lo = 0x0000000000000000}},
        M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(C[i].v[0], ==, R[i].v[0]);
        assert_ullong(C[i].v[1], ==, R[i].v[1]);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_qgemm_3x2x1_0(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float128_t alpha = { SB_REAL128L_ONE, SB_REAL128U_ONE };
    const float128_t beta = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    const uint64_t M = 3;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 1;     // Number of columns in B

    // 3x2 matrix A (row-major order)
    float128_t* A = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000}},
        M*N);

    // 2x1 matrix B (row-major order)
    float128_t* B = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000}},
        N*P);

    // 3x1 matrix C (row-major order)
    float128_t* C = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call qgemm
    qgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float128_t* R = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002600000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003100000000000, .lo = 0x0000000000000000}},
        M*P);

    for (uint64_t i = 0; i < 3; i++) {
        assert_ullong(C[i].v[0], ==, R[i].v[0]);
        assert_ullong(C[i].v[1], ==, R[i].v[1]);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_qgemm_5x4x3(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float128_t alpha = { 0x0000000000000000, 0x4000400000000000 };
    const float128_t beta = { SB_REAL128L_NEGONE, SB_REAL128U_NEGONE };

    const uint64_t M = 5;     // Number of rows in A
    const uint64_t N = 4;     // Number of columns in A
    const uint64_t P = 3;     // Number of columns in B

    // 5x4 matrix A (row-major order)
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

    // 4x3 matrix B (row-major order)
    float128_t* B = qvec((float128_pair_t[]){
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002000000000000, .lo = 0x0000000000000000},

            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001c00000000000, .lo = 0x0000000000000000},

            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},

            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000}},
        N*P);

    // 5x3 matrix C (row-major order)
    float128_t* C = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},

            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000},

            {.hi = 0xbfff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000800000000000, .lo = 0x0000000000000000},

            {.hi = 0x1000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x1400000000000000, .lo = 0x0000000000000000},
            {.hi = 0x1800800000000000, .lo = 0x0000000000000000},

            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000}},
        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call qgemm
    qgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float128_t* R = qvec((float128_pair_t[]){
            {.hi = 0x4002700000000000, .lo = 0x0000000000000000},
            {.hi = 0x3ffe000000000000, .lo = 0x0000000000000000},
            {.hi = 0x40048c0000000000, .lo = 0x0000000000000000},

            {.hi = 0x40056c0000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005960000000000, .lo = 0x0000000000000000},
            
            {.hi = 0x4005120000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005520000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005920000000000, .lo = 0x0000000000000000},
            
            {.hi = 0x4005220000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003e00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4005180000000000, .lo = 0x0000000000000000},
            
            {.hi = 0x4006030000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004f80000000000, .lo = 0x0000000000000000},
            {.hi = 0x40063e0000000000, .lo = 0x0000000000000000}},
        M*P);

    for (uint64_t i = 0; i < M*P; i++) {
        assert_ullong(C[i].v[0], ==, R[i].v[0]);
        assert_ullong(C[i].v[1], ==, R[i].v[1]);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}
