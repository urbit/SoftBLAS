#include "test.h"
#include <stdio.h>

MunitResult test_sgemm_0_row(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float32_t alpha = { SB_REAL32_ONE };  // Scalar alpha
    const float32_t beta = { SB_REAL32_ZERO };  // Scalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float32_t* A = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f}, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float32_t* B = svec((float[]){5.0f, 6.0f, 7.0f, 8.0f}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float32_t* C = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f}, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call sgemm
    sgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float32_t* R = svec((float[]){19.0f, 22.0f, 43.0f, 50.0f}, M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_sgemm_0_col(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'T';  // Transpose
    const char transB = 'T';  // Transpose
    const float32_t alpha = { SB_REAL32_ONE };  // Scalar alpha
    const float32_t beta = { SB_REAL32_ZERO };  // Scalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float32_t* A = svec((float[]){1.0f, 3.0f, 2.0f, 4.0f}, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float32_t* B = svec((float[]){5.0f, 7.0f, 6.0f, 8.0f}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float32_t* C = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f}, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call sgemm
    sgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float32_t* R = svec((float[]){19.0f, 22.0f, 43.0f, 50.0f}, M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_sgemm_3x2x1_0(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'T';  // No transpose
    const char transB = 'T';  // No transpose
    const float32_t alpha = {*(uint32_t*)&(float){2.5f}};
    const float32_t beta = {*(uint32_t*)&(float){0.0f}};

    const uint64_t M = 3;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 1;     // Number of columns in B

    // 3x2 matrix A (row-major order)
    float32_t* A = svec((float[]){ 1.0f,  2.0f,
                                   3.0f,  4.0f,
                                   5.0f,  6.0f},
                        M*N);

    // 2x1 matrix B (row-major order)
    float32_t* B = svec((float[]){ 1.0f,
                                   2.0f},
                        N*P);

    // 3x1 matrix C (row-major order)
    float32_t* C = svec((float[]){ 0.0f,
                                   0.0f,
                                   0.0f},
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call sgemm
    sgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float32_t* R = svec((float[]){ 12.5f,
                                   27.5f,
                                   42.5f},
                        M*P);

    for (uint64_t i = 0; i < M*P; i++) {
        printf("C[%d] = %x\tR[%d] = %x\n", i, C[i].v, i, R[i].v);
    }

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_sgemm_5x4x3(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';  // No transpose
    const char transB = 'N';  // No transpose
    const float32_t alpha = {*(uint32_t*)&(float){2.5f}};
    const float32_t beta = {*(uint32_t*)&(float){0.0f}};

    const uint64_t M = 5;     // Number of rows in A
    const uint64_t N = 4;     // Number of columns in A
    const uint64_t P = 3;     // Number of columns in B

    // 5x4 matrix A (row-major order)
    float32_t* A = svec((float[]){ 2.0f, -1.0f,  3.0f,  0.0f,
                                   4.0f,  1.0f,  0.0f,  2.0f,
                                  -1.0f,  3.0f,  4.0f,  5.0f,
                                   1.0f,  2.0f,  0.0f,  3.0f,
                                   2.0f,  4.0f,  2.0f,  5.0f},
                        M*N);

    // 4x3 matrix B (row-major order)
    float32_t* B = svec((float[]){ 6.0f, -2.0f,  8.0f,
                                   4.0f,  1.0f,  7.0f,
                                  -1.0f,  2.0f,  4.0f,
                                   5.0f,  4.0f,  2.0f},
                        N*P);

    // 2x2 matrix C (row-major order)
    float32_t* C = svec((float[]){ 1.0f,  2.0f,  3.0f,
                                   4.0f,  5.0f,  6.0f,
                                  -1.0f, -2.0f, -3.0f,
                                  -4.0f, -5.0f, -6.0f,
                                  -2.0f, -3.0f, -4.0f},
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call sgemm
    sgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float32_t* R = svec((float[]){11.5f,  0.5f,  49.5f,
                                  91.0f, -2.5f, 101.5f,
                                  68.5f, 84.5f, 100.5f,
                                  76.5f, 35.0f,  76.0f,
                                 129.0f, 63.0f, 159.0f},
                        M*P);

    for (uint64_t i = 0; i < M*P; i++) {
        printf("C[%d] = %x\tR[%d] = %x\n", i, C[i].v, i, R[i].v);
    }

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}
