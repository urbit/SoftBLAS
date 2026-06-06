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
