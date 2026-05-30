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
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY, 'n');

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
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY, 'n');

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
    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY, 'n');

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

    // 4x5 matrix A (row-major order); only vector Y is strided (incY=2),
    // so A is stored densely with no padding rows.
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
    float32_t* SY = svec((float[]){-3.0f, 0.0f, 4.0f, 0.0f, -5.0f, 0.0f, 6.0f}, 7);

    // Stride for vector X and Y
    const int64_t incX = 1;
    const int64_t incY = 2;

    // Leading dimension of A
    const uint64_t lda = N;

    // Call sgemv
    sgemv(Layout, Trans, 4, 5, alpha, A, lda, SX, incX, beta, SY, incY, 'n');

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

//  Regression for the gemv stride bug: X is strided (incX=2) while the matrix
//  A is dense. The old code folded incX into the A index and gave wrong results
//  (and could read out of bounds) for any incX != 1.
MunitResult test_sgemv_incx(const MunitParameter params[],
                            void *user_data) {
    const char Layout = 'R';  // Row-major
    const char Trans = 'N';   // No transpose
    const float32_t alpha = { SB_REAL32_ONE };
    const float32_t beta = { SB_REAL32_ZERO };

    const uint64_t M = 4;
    const uint64_t N = 5;
    float32_t* A = svec((float[]){ 2.0f, -1.0f,  3.0f,  0.0f,  4.0f,
                                   1.0f,  0.0f,  2.0f, -1.0f,  3.0f,
                                   4.0f,  5.0f,  1.0f,  2.0f,  0.0f,
                                   3.0f,  2.0f,  4.0f,  2.0f,  5.0f},
                        M*N);

    // X = [1 2 3 4 5] stored with stride 2 (padding entries are never read).
    float32_t* SX = svec((float[]){1.0f, 99.0f, 2.0f, 99.0f, 3.0f,
                                   99.0f, 4.0f, 99.0f, 5.0f}, 9);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f}, 4);

    const int64_t incX = 2;
    const int64_t incY = 1;
    const uint64_t lda = N;

    sgemv(Layout, Trans, M, N, alpha, A, lda, SX, incX, beta, SY, incY, 'n');

    float32_t* RY = svec((float[]){29.0f, 18.0f, 25.0f, 52.0f}, 4);

    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(A); free(SX); free(SY); free(RY);
    return MUNIT_OK;
}

//  Regression: an invalid Layout must return without modifying Y or killing
//  the process (these routines previously called exit(-1) on bad arguments).
MunitResult test_sgemv_badlayout(const MunitParameter params[],
                                 void *user_data) {
    const float32_t alpha = { SB_REAL32_ONE };
    const float32_t beta = { SB_REAL32_ONE };
    float32_t* A = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f}, 4);
    float32_t* X = svec((float[]){1.0f, 2.0f}, 2);
    float32_t* Y = svec((float[]){5.0f, 6.0f}, 2);

    //  'X' is not a valid Layout; the routine must return immediately.
    sgemv('X', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');

    //  Y must be untouched.
    float32_t* RY = svec((float[]){5.0f, 6.0f}, 2);
    for (uint64_t i = 0; i < 2; i++) {
        assert_ulong(Y[i].v, ==, RY[i].v);
    }

    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}

//  Transpose path (Trans='T'): y = A^T x. A is 2x3 row-major [[1,2,3],[4,5,6]],
//  x=[1,1] (len M=2), so y = [1+4, 2+5, 3+6] = [5,7,9] (len N=3). No gemv test
//  exercised the transpose branch before.
MunitResult test_sgemv_trans(const MunitParameter params[], void *user_data) {
    const float32_t alpha = { SB_REAL32_ONE }, beta = { SB_REAL32_ZERO };
    float32_t* A = svec((float[]){1.0f,2.0f,3.0f, 4.0f,5.0f,6.0f}, 6);
    float32_t* X = svec((float[]){1.0f, 1.0f}, 2);
    float32_t* Y = svec((float[]){0.0f, 0.0f, 0.0f}, 3);
    sgemv('R', 'T', 2, 3, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float32_t* RY = svec((float[]){5.0f, 7.0f, 9.0f}, 3);
    for (uint64_t i = 0; i < 3; i++) assert_ulong(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}

//  Padded leading dimension: A is logically 2x2 [[1,2],[3,4]] but stored with
//  lda=3 (a trailing pad column). y = A*[1,1] = [3,7]; the pad must be skipped.
MunitResult test_sgemv_padlda(const MunitParameter params[], void *user_data) {
    const float32_t alpha = { SB_REAL32_ONE }, beta = { SB_REAL32_ZERO };
    float32_t* A = svec((float[]){1.0f,2.0f,99.0f, 3.0f,4.0f,99.0f}, 6);
    float32_t* X = svec((float[]){1.0f, 1.0f}, 2);
    float32_t* Y = svec((float[]){0.0f, 0.0f}, 2);
    sgemv('R', 'N', 2, 2, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float32_t* RY = svec((float[]){3.0f, 7.0f}, 2);
    for (uint64_t i = 0; i < 2; i++) assert_ulong(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
