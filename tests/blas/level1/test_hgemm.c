#include "test.h"
#include <stdio.h>

MunitResult test_hgemm_0_row(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';   // No transpose
    const char transB = 'N';   // No transpose
    const float16_t alpha = { SB_REAL16_ONE };  // Hcalar alpha
    const float16_t beta = { SB_REAL16_ZERO };  // Hcalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float16_t* A = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 }, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float16_t* B = hvec((uint16_t[]){0x4500 , 0x4600, 0x4700, 0x4800}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float16_t* C = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 }, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call hgemm
    hgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float16_t* R = hvec((uint16_t[]){0x4cc0, 0x4d80 , 0x5160 , 0x5240 }, M*P);

    for (uint64_t i = 0; i < M*P; i++) {
        assert_ushort(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_hgemm_0_col(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'T';  // Transpose
    const char transB = 'T';  // Transpose
    const float16_t alpha = { SB_REAL16_ONE };  // Hcalar alpha
    const float16_t beta = { SB_REAL16_ZERO };  // Hcalar beta

    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 2;     // Number of columns in B

    // 2x2 matrix A (row-major order)
    //   [1.0 2.0]
    //   [3.0 4.0]
    float16_t* A = hvec((uint16_t[]){0x3c00 , 0x4200 , 0x4000 , 0x4400 }, M*N);

    // 2x2 matrix B (row-major order)
    //   [5.0 6.0]
    //   [7.0 8.0]
    float16_t* B = hvec((uint16_t[]){0x4500 , 0x4700, 0x4600, 0x4800}, N*P);

    // 2x2 matrix C (row-major order)
    //   [0.0 0.0]
    //   [0.0 0.0]
    float16_t* C = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 }, M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call hgemm
    hgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float16_t* R = hvec((uint16_t[]){0x4cc0, 0x4d80 , 0x5160 , 0x5240 }, M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ushort(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_hgemm_3x2x1_0(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';  // No transpose
    const char transB = 'N';  // No transpose
    const float16_t alpha = {*(uint16_t*)&(uint16_t){0x3c00 }};
    const float16_t beta = {*(uint16_t*)&(uint16_t){0x0 }};

    const uint64_t M = 3;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    const uint64_t P = 1;     // Number of columns in B

    // 3x2 matrix A (row-major order)
    float16_t* A = hvec((uint16_t[]){ 0x3c00 ,  0x4000 ,
                                   0x4200 ,  0x4400 ,
                                   0x4500 ,  0x4600},
                        M*N);

    // 2x1 matrix B (row-major order)
    float16_t* B = hvec((uint16_t[]){ 0x3c00 ,
                                   0x4000 },
                        N*P);

    // 3x1 matrix C (row-major order)
    float16_t* C = hvec((uint16_t[]){ 0x0 ,
                                   0x0 ,
                                   0x0 },
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call hgemm
    hgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float16_t* R = hvec((uint16_t[]){ 0x4500 ,
                                   0x4980 ,
                                   0x4c40},
                        M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ushort(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}

MunitResult test_hgemm_5x4x3(const MunitParameter params[],
                             void *user_data) {
    const char transA = 'N';  // No transpose
    const char transB = 'N';  // No transpose
    const float16_t alpha = {*(uint16_t*)&(uint16_t){0x4100}};
    const float16_t beta = {*(uint16_t*)&(uint16_t){0xbc00 }};

    const uint64_t M = 5;     // Number of rows in A
    const uint64_t N = 4;     // Number of columns in A
    const uint64_t P = 3;     // Number of columns in B

    // 5x4 matrix A (row-major order)
    float16_t* A = hvec((uint16_t[]){ 0x4000 , 0xbc00 ,  0x4200 ,  0x0 ,
                                   0x4400 ,  0x3c00 ,  0x0 ,  0x4000 ,
                                  0xbc00 ,  0x4200 ,  0x4400 ,  0x4500 ,
                                   0x3c00 ,  0x4000 ,  0x0 ,  0x4200 ,
                                   0x4000 ,  0x4400 ,  0x4000 ,  0x4500 },
                        M*N);

    // 4x3 matrix B (row-major order)
    float16_t* B = hvec((uint16_t[]){ 0x4600, 0xc000 ,  0x4800,
                                   0x4400 ,  0x3c00 ,  0x4700,
                                  0xbc00 ,  0x4000 ,  0x4400 ,
                                   0x4500 ,  0x4400 ,  0x4000 },
                        N*P);

    // 2x2 matrix C (row-major order)
    float16_t* C = hvec((uint16_t[]){ 0x3c00 ,  0x4000 ,  0x4200 ,
                                   0x4400 ,  0x4500 ,  0x4600,
                                  0xbc00 , 0xc000 , 0xc200 ,
                                  0xc400 , 0xc500 , 0xc600,
                                  0xc000 , 0xc200 , 0xc400 },
                        M*P);

    // Leading dimensions
    const uint64_t lda = N;
    const uint64_t ldb = P;
    const uint64_t ldc = P;

    // Call hgemm
    hgemm(transA, transB, M, N, P, alpha, A, lda, B, ldb, beta, C, ldc);

    float16_t* R = hvec((uint16_t[]){0x49c0, 0x3800, 0x5230,
                                     0x55b0, 0xc100, 0x5658,
                                     0x5448, 0x5548, 0x5648,
                                     0x54c8, 0x5060, 0x54c0,
                                     0x580c, 0x53e0, 0x58f8},
                        M*P);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ushort(C[i].v, ==, R[i].v);
    }

    free(A);
    free(B);
    free(C);
    free(R);

    return MUNIT_OK;
}