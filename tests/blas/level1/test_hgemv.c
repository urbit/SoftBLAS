#include "test.h"

MunitResult test_hgemv_0_row(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float16_t alpha = { SB_REAL16_ONE };  // Hcalar alpha
    const float16_t beta = { SB_REAL16_ZERO };  // Hcalar beta

    // 2x2 matrix A (row-major order)
    //   [1.0 3.0]
    //   [2.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float16_t* A = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 }, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 }, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float16_t* HY = hvec((uint16_t[]){0x4200 , 0x4400 }, 2);

    // Htride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call hgemv
    hgemv(Layout, Trans, M, N, alpha, A, lda, HX, incX, beta, HY, incY);

    float16_t* RY = hvec((uint16_t[]){0x4500 , 0x4980 }, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(A);
    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_hgemv_0_col(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'C';  // Column-major order
    const char Trans = 'N';   // No transpose
    const float16_t alpha = { SB_REAL16_ONE };  // Hcalar alpha
    const float16_t beta = { SB_REAL16_ZERO };  // Hcalar beta

    // 2x2 matrix A (column-major order)
    //   [1.0 3.0]
    //   [2.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float16_t* A = hvec((uint16_t[]){0x3c00 , 0x4200 , 0x4000 , 0x4400 }, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 }, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float16_t* HY = hvec((uint16_t[]){0x4200 , 0x4400 }, 2);

    // Htride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = M;

    // Call hgemv
    hgemv(Layout, Trans, M, N, alpha, A, lda, HX, incX, beta, HY, incY);

    float16_t* RY = hvec((uint16_t[]){0x4500 , 0x4980 }, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(A);
    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_hgemv_12345(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float16_t alpha = {*(uint16_t*)&(uint16_t){0x4100}};
    const float16_t beta = {*(uint16_t*)&(uint16_t){0xbc00 }};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 4;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float16_t* A = hvec((uint16_t[]){ 0x4000 , 0xbc00 ,  0x4200 ,  0x0 ,  0x4400 ,
                                   0x3c00 ,  0x0 ,  0x4000 , 0xbc00 ,  0x4200 ,
                                   0x4400 ,  0x4500 ,  0x3c00 ,  0x4000 ,  0x0 ,
                                   0x4200 ,  0x4000 ,  0x4400 ,  0x4000 ,  0x4500 },
                        M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float16_t* HY = hvec((uint16_t[]){0xc200 , 0x4400 , 0xc500 , 0x4600}, 4);

    // Htride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call hgemv
    hgemv(Layout, Trans, M, N, alpha, A, lda, HX, incX, beta, HY, incY);

    float16_t* RY = hvec((uint16_t[]){0x54b8, 0x5120 , 0x5438, 0x57c0 }, 4);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(A);
    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_hgemv_stride(const MunitParameter params[],
                              void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float16_t alpha = {*(uint16_t*)&(uint16_t){0x4100}};
    const float16_t beta = {*(uint16_t*)&(uint16_t){0xbc00 }};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 7;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float16_t* A = hvec((uint16_t[]){ 0x4000 , 0xbc00 ,  0x4200 ,  0x0 ,  0x4400 ,
                                   0x0 ,  0x0 ,  0x0 ,  0x0 ,  0x0 ,
                                   0x3c00 ,  0x0 ,  0x4000 , 0xbc00 ,  0x4200 ,
                                   0x0 ,  0x0 ,  0x0 ,  0x0 ,  0x0 ,
                                   0x4400 ,  0x4500 ,  0x3c00 ,  0x4000 ,  0x0 ,
                                   0x0 ,  0x0 ,  0x0 ,  0x0 ,  0x0 ,
                                   0x4200 ,  0x4000 ,  0x4400 ,  0x4000 ,  0x4500 },
                       M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float16_t* HY = hvec((uint16_t[]){0xc200 , 0x0 , 0x4400 , 0x0 , 0xc500 , 0x0 , 0x4600}, 7);

    // Htride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 2;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call hgemv
    hgemv(Layout, Trans, 4, 5, alpha, A, lda, HX, incX, beta, HY, incY);

    float16_t* RY = hvec((uint16_t[]){0x54b8, 0x0 , 0x5120 , 0x0 , 0x5438, 0x0 , 0x57c0 }, 7);

    for (uint64_t i = 0; i < 7; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(A);
    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}
