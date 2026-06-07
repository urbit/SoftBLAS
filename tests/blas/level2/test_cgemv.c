#include "test.h"

//  Complex GEMV (half 'i' / single 'c' / double 'z' / quad 'v').
//  y <- alpha*op(A)*x + beta*y, op = N/T/C. Layout 'R'(row) / 'C'(col).
//  Small integer values keep every result exact and bit-checkable.

//  Quad helpers (high word holds the value; low word stays 0).
#define QH(hi) {{ 0x0ull, (hi) }}
#define Q_ONE   0x3fff000000000000ull
#define Q_TWO   0x4000000000000000ull
#define Q_THREE 0x4000800000000000ull

//  Single, row-major, no transpose. A=[[1+1i, 2],[0, 1]], x=[1, 1+1i] ->
//  y = [ (1+1i)*1 + 2*(1+1i), 1*(1+1i) ] = [3+3i, 1+1i].
MunitResult test_cgemv_basic(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,0, 0,0, 1,0}, 4);
    complex32_t* X = cvec((float[]){1,0, 1,1}, 2);
    complex32_t* Y = cvec((float[]){0,0, 0,0}, 2);
    cgemv('R', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ulong(Y[0].real.v, ==, 0x40400000u); assert_ulong(Y[0].imag.v, ==, 0x40400000u);  // 3+3i
    assert_ulong(Y[1].real.v, ==, 0x3f800000u); assert_ulong(Y[1].imag.v, ==, 0x3f800000u);  // 1+1i
    free(A); free(X); free(Y);
    return MUNIT_OK;
}

//  Same logical operation, column-major storage (A stored by columns, lda=M).
MunitResult test_cgemv_colmajor(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 0,0, 2,0, 1,0}, 4);   // columns: [1+1i,0],[2,1]
    complex32_t* X = cvec((float[]){1,0, 1,1}, 2);
    complex32_t* Y = cvec((float[]){0,0, 0,0}, 2);
    cgemv('C', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ulong(Y[0].real.v, ==, 0x40400000u); assert_ulong(Y[0].imag.v, ==, 0x40400000u);  // 3+3i
    assert_ulong(Y[1].real.v, ==, 0x3f800000u); assert_ulong(Y[1].imag.v, ==, 0x3f800000u);  // 1+1i
    free(A); free(X); free(Y);
    return MUNIT_OK;
}

//  Conjugate transpose: y = conj(A^T)*x. A=[[1+1i, 2+1i],[3-1i, 0+2i]], x=[1,1]
//  -> conj(A^T)=[[1-1i, 3+1i],[2-1i, 0-2i]], y = [4+0i, 2-3i].
MunitResult test_cgemv_conjtrans(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,1, 3,-1, 0,2}, 4);
    complex32_t* X = cvec((float[]){1,0, 1,0}, 2);
    complex32_t* Y = cvec((float[]){0,0, 0,0}, 2);
    cgemv('R', 'C', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ulong(Y[0].real.v, ==, 0x40800000u); assert_ulong(Y[0].imag.v, ==, 0x00000000u);  // 4+0i
    assert_ulong(Y[1].real.v, ==, 0x40000000u); assert_ulong(Y[1].imag.v, ==, 0xc0400000u);  // 2-3i
    free(A); free(X); free(Y);
    return MUNIT_OK;
}

//  Double, row-major, no transpose (same matrices as test_cgemv_basic).
MunitResult test_zgemv_basic(const MunitParameter params[], void* u) {
    const complex64_t alpha = {{ 0x3ff0000000000000ull }, { 0 }};  // 1
    const complex64_t beta  = {{ 0 }, { 0 }};                       // 0
    complex64_t A[4] = {{{ 0x3ff0000000000000ull }, { 0x3ff0000000000000ull }},  // 1+1i
                        {{ 0x4000000000000000ull }, { 0 }},                       // 2
                        {{ 0 }, { 0 }},                                           // 0
                        {{ 0x3ff0000000000000ull }, { 0 }}};                      // 1
    complex64_t X[2] = {{{ 0x3ff0000000000000ull }, { 0 }},                       // 1
                        {{ 0x3ff0000000000000ull }, { 0x3ff0000000000000ull }}};  // 1+1i
    complex64_t Y[2] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    zgemv('R', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ullong(Y[0].real.v, ==, 0x4008000000000000ull); assert_ullong(Y[0].imag.v, ==, 0x4008000000000000ull);  // 3+3i
    assert_ullong(Y[1].real.v, ==, 0x3ff0000000000000ull); assert_ullong(Y[1].imag.v, ==, 0x3ff0000000000000ull);  // 1+1i
    return MUNIT_OK;
}

//  Half, row-major, no transpose (same matrices as test_cgemv_basic).
MunitResult test_igemv_basic(const MunitParameter params[], void* u) {
    const complex16_t alpha = SB_COMPLEX16_ONE, beta = SB_COMPLEX16_ZERO;
    complex16_t A[4] = {{{ 0x3c00 }, { 0x3c00 }}, {{ 0x4000 }, { 0 }},
                        {{ 0 }, { 0 }},           {{ 0x3c00 }, { 0 }}};   // [[1+1i,2],[0,1]]
    complex16_t X[2] = {{{ 0x3c00 }, { 0 }}, {{ 0x3c00 }, { 0x3c00 }}};   // [1, 1+1i]
    complex16_t Y[2] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    igemv('R', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ushort(Y[0].real.v, ==, 0x4200); assert_ushort(Y[0].imag.v, ==, 0x4200);  // 3+3i
    assert_ushort(Y[1].real.v, ==, 0x3c00); assert_ushort(Y[1].imag.v, ==, 0x3c00);  // 1+1i
    return MUNIT_OK;
}

//  Quad, row-major, no transpose (same matrices as test_cgemv_basic).
MunitResult test_vgemv_basic(const MunitParameter params[], void* u) {
    const complex128_t alpha = { QH(Q_ONE), QH(0x0ull) };  // 1
    const complex128_t beta  = { QH(0x0ull), QH(0x0ull) };  // 0
    complex128_t A[4] = {{ QH(Q_ONE), QH(Q_ONE) }, { QH(Q_TWO), QH(0x0ull) },
                         { QH(0x0ull), QH(0x0ull) }, { QH(Q_ONE), QH(0x0ull) }};  // [[1+1i,2],[0,1]]
    complex128_t X[2] = {{ QH(Q_ONE), QH(0x0ull) }, { QH(Q_ONE), QH(Q_ONE) }};    // [1, 1+1i]
    complex128_t Y[2] = {{ QH(0x0ull), QH(0x0ull) }, { QH(0x0ull), QH(0x0ull) }};
    vgemv('R', 'N', 2, 2, alpha, A, 2, X, 1, beta, Y, 1, 'n');
    assert_ullong(Y[0].real.v[1], ==, Q_THREE); assert_ullong(Y[0].imag.v[1], ==, Q_THREE);  // 3+3i
    assert_ullong(Y[1].real.v[1], ==, Q_ONE);   assert_ullong(Y[1].imag.v[1], ==, Q_ONE);    // 1+1i
    return MUNIT_OK;
}

//  3x3 single (row-major, N). A=A3, x=[1, 1+1i, 2-1i] -> y=[4+5i, 6+1i, 6+3i].
MunitResult test_cgemv_3x3(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,0, 0,1, 0,2, 1,1, 3,0, 1,0, 0,1, 2,2}, 9);
    complex32_t* X = cvec((float[]){1,0, 1,1, 2,-1}, 3);
    complex32_t* Y = cvec((float[]){0,0, 0,0, 0,0}, 3);
    cgemv('R', 'N', 3, 3, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    complex32_t* R = cvec((float[]){4,5, 6,1, 6,3}, 3);
    for (uint64_t i = 0; i < 3; i++) {
        assert_ulong(Y[i].real.v, ==, R[i].real.v);
        assert_ulong(Y[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(X); free(Y); free(R);
    return MUNIT_OK;
}

//  4x4 double (row-major, N). x=[1, i, 2, 1-i] -> y=[3+4i, 5-1i, 6+4i, 5-4i].
MunitResult test_zgemv_4x4(const MunitParameter params[], void* u) {
    const complex64_t alpha = {{ 0x3ff0000000000000ull }, { 0 }}, beta = {{ 0 }, { 0 }};
    complex64_t* A = zvec((double[]){1,1, 2,-1, 0,1, 1,0, 0,2, 1,1, 3,0, 2,-2,
                                     1,0, 0,1, 2,2, 1,1, 3,-1, 1,0, 0,-1, 2,0}, 16);
    complex64_t* X = zvec((double[]){1,0, 0,1, 2,0, 1,-1}, 4);
    complex64_t* Y = zvec((double[]){0,0, 0,0, 0,0, 0,0}, 4);
    zgemv('R', 'N', 4, 4, alpha, A, 4, X, 1, beta, Y, 1, 'n');
    complex64_t* R = zvec((double[]){3,4, 5,-1, 6,4, 5,-4}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ullong(Y[i].real.v, ==, R[i].real.v);
        assert_ullong(Y[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(X); free(Y); free(R);
    return MUNIT_OK;
}

//  25x25 single with vector stride incX=incY=5. A and x are generated from a
//  deterministic (non-periodic, small-integer) formula so every value is exact;
//  the expected y is computed independently (see the Python oracle). The skipped
//  X slots hold a 99+99i sentinel, so a stride-indexing bug reads garbage and fails.
MunitResult test_cgemv_25x25_stride5(const MunitParameter params[], void* u) {
    const uint64_t N = 25;
    const int64_t inc = 5;
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;

    complex32_t* A = malloc(N * N * sizeof(complex32_t));
    for (uint64_t i = 0; i < N; i++) {
        for (uint64_t j = 0; j < N; j++) {
            float re = (float)((int)((i * 7 + j * j) % 11) - 5);
            float im = (float)((int)((i * 3 + j * 5) % 13) - 6);
            A[i * N + j].real.v = *(uint32_t*)&re;
            A[i * N + j].imag.v = *(uint32_t*)&im;
        }
    }
    complex32_t* X = malloc(N * (uint64_t)inc * sizeof(complex32_t));
    for (uint64_t k = 0; k < N * (uint64_t)inc; k++) {
        float s = 99.0f;
        X[k].real.v = *(uint32_t*)&s;
        X[k].imag.v = *(uint32_t*)&s;
    }
    for (uint64_t j = 0; j < N; j++) {
        float re = (float)((int)(j % 3) - 1);
        float im = (float)((int)((2 * j + 1) % 3) - 1);
        X[j * (uint64_t)inc].real.v = *(uint32_t*)&re;
        X[j * (uint64_t)inc].imag.v = *(uint32_t*)&im;
    }
    complex32_t* Y = malloc(N * (uint64_t)inc * sizeof(complex32_t));
    for (uint64_t k = 0; k < N * (uint64_t)inc; k++) { Y[k].real.v = 0; Y[k].imag.v = 0; }

    cgemv('R', 'N', N, N, alpha, A, N, X, inc, beta, Y, inc, 'n');

    complex32_t* R = cvec((float[]){
        -8,23, 22,7, 4,-20, -7,1, 10,22, -8,-3, 9,-6, -11,-22, -9,12, -3,11,
        3,-18, 5,3, 9,0, -22,12, 19,7, 23,2, -21,-10, -4,0, -11,-3, 17,5,
        -3,-22, -23,1, 5,22, 22,-7, -9,-8}, 25);
    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(Y[i * (uint64_t)inc].real.v, ==, R[i].real.v);
        assert_ulong(Y[i * (uint64_t)inc].imag.v, ==, R[i].imag.v);
    }
    free(A); free(X); free(Y); free(R);
    return MUNIT_OK;
}

//  NaN canonicalization through cgemv: a NaN in A taints Y[0]; nan_unify_c
//  forces both components to SINGNAN.
MunitResult test_cgemv_nan(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t A[1] = {{{ 0x7fc00001u }, { 0 }}};   // NaN(non-canon) + 0i
    complex32_t X[1] = {{{ 0x3f800000u }, { 0 }}};   // 1 + 0i
    complex32_t Y[1] = {{{ 0 }, { 0 }}};
    cgemv('R', 'N', 1, 1, alpha, A, 1, X, 1, beta, Y, 1, 'n');
    assert_ulong(Y[0].real.v, ==, (uint32_t)SINGNAN);
    assert_ulong(Y[0].imag.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
