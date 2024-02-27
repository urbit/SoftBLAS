#include "test.h"

MunitResult test_sgemv_0_row(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float32_t alpha = { SB_REAL32_ONE };  // Scalar alpha
    const float32_t beta = { SB_REAL32_ZERO };  // Scalar beta

    // 2x2 matrix A (row-major order)
    //   [1.0 3.0]
    //   [2.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float32_t* A = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f}, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float32_t* SX = svec((float[]){1.0f, 2.0f}, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float32_t* SY = svec((float[]){3.0f, 4.0f}, 2);

    // Stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call sgemv
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY);

    float32_t* RY = svec((float[]){5.0f, 11.0f}, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(A);
    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_sgemv_0_col(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'C';  // Column-major order
    const char Trans = 'N';   // No transpose
    const float32_t alpha = { SB_REAL32_ONE };  // Scalar alpha
    const float32_t beta = { SB_REAL32_ZERO };  // Scalar beta

    // 2x2 matrix A (column-major order)
    //   [1.0 3.0]
    //   [2.0 4.0]
    const uint64_t M = 2;     // Number of rows in A
    const uint64_t N = 2;     // Number of columns in A
    float32_t* A = svec((float[]){1.0f, 3.0f, 2.0f, 4.0f}, 4);

    // Vector X (column vector)
    //   [1.0]
    //   [2.0]
    float32_t* SX = svec((float[]){1.0f, 2.0f}, 2);

    // Vector Y (column vector)
    //   [3.0]
    //   [4.0]
    float32_t* SY = svec((float[]){3.0f, 4.0f}, 2);

    // Stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = M;

    // Call sgemv
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY);

    float32_t* RY = svec((float[]){5.0f, 11.0f}, 2);

    for (uint64_t i = 0; i < 2; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(A);
    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_sgemv_12345(const MunitParameter params[],
                             void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float32_t alpha = {*(uint32_t*)&(float){2.5f}};
    const float32_t beta = {*(uint32_t*)&(float){-1.0f}};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 4;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float32_t* A = svec((float[]){ 2.0f, -1.0f,  3.0f,  0.0f,  4.0f,
                                   1.0f,  0.0f,  2.0f, -1.0f,  3.0f,
                                   4.0f,  5.0f,  1.0f,  2.0f,  0.0f,
                                   3.0f,  2.0f,  4.0f,  2.0f,  5.0f},
                        M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float32_t* SY = svec((float[]){-3.0f, 4.0f, -5.0f, 6.0f}, 4);

    // Stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 1;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call sgemv
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY);

    float32_t* RY = svec((float[]){75.5f, 41.0f, 67.5f, 124.0f}, 4);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(A);
    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_sgemv_stride(const MunitParameter params[],
                              void *user_data) {
    const char Layout = 'R';  // Row-major order
    const char Trans = 'N';   // No transpose
    const float32_t alpha = {*(uint32_t*)&(float){2.5f}};
    const float32_t beta = {*(uint32_t*)&(float){-1.0f}};

    // 5x4 matrix A (row-major order)
    const uint64_t M = 7;     // Number of rows in A
    const uint64_t N = 5;     // Number of columns in A
    float32_t* A = svec((float[]){ 2.0f, -1.0f,  3.0f,  0.0f,  4.0f,
                                   0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
                                   1.0f,  0.0f,  2.0f, -1.0f,  3.0f,
                                   0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
                                   4.0f,  5.0f,  1.0f,  2.0f,  0.0f,
                                   0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
                                   3.0f,  2.0f,  4.0f,  2.0f,  5.0f},
                       M*N);

    // Vector X (column vector)
    //   [1.0 2.0 3.0 4.0 5.0]'
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);

    // Vector Y (column vector)
    //   [-3.0 4.0 -5.0 6.0]'
    float32_t* SY = svec((float[]){-3.0f, 0.0f, 4.0f, 0.0f, -5.0f, 0.0f, 6.0f}, 7);

    // Stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 2;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call sgemv
    sgemv(Layout, Trans, 4, 5, alpha, A, lda, SX, incX, beta, SY, incY);

    float32_t* RY = svec((float[]){75.5f, 0.0f, 41.0f, 0.0f, 67.5f, 0.0f, 124.0f}, 7);

    for (uint64_t i = 0; i < 7; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(A);
    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}
