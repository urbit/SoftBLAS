#include "test.h"

//  Complex GEMM (single 'c' / double 'z'). C <- alpha*op(A)*op(B) + beta*C,
//  op = N/T/C. Row-major, leading dims lda/ldb/ldc, mirroring sgemm. All test
//  values are small integers so the float results are exact and bit-checkable.

//  Canonical vector from the brief:
//  A=[[1+1i, 2],[0, 1]], B=[[1, 0],[0, 1+1i]] -> C=[[1+1i, 2+2i],[0, 1+1i]].
MunitResult test_cgemm_basic(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE;
    const complex32_t beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,0, 0,0, 1,0}, 4);
    complex32_t* B = cvec((float[]){1,0, 0,0, 0,0, 1,1}, 4);
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0}, 4);
    cgemm('N', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    complex32_t* R = cvec((float[]){1,1, 2,2, 0,0, 1,1}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  transA='C' is the conjugate transpose. A stored [[1+1i, 2+1i],[3-1i, 0+2i]],
//  so op(A)=conj(A^T)=[[1-1i, 3+1i],[2-1i, 0-2i]]; times identity B (transB='N')
//  reproduces op(A).
MunitResult test_cgemm_conjtrans(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE;
    const complex32_t beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,1, 3,-1, 0,2}, 4);
    complex32_t* B = cvec((float[]){1,0, 0,0, 0,0, 1,0}, 4);   // identity
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0}, 4);
    cgemm('C', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    complex32_t* R = cvec((float[]){1,-1, 3,1, 2,-1, 0,-2}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  General alpha/beta (1x1): alpha*(A*B) + beta*C with alpha=1+1i, beta=2,
//  A=2, B=3+4i, C=1+1i -> (1+1i)(6+8i) + 2(1+1i) = (-2+14i)+(2+2i) = 0 + 16i.
MunitResult test_cgemm_alphabeta(const MunitParameter params[], void* u) {
    const complex32_t alpha = {{ 0x3f800000 }, { 0x3f800000 }};  // 1 + 1i
    const complex32_t beta  = {{ 0x40000000 }, { 0x00000000 }};  // 2 + 0i
    complex32_t* A = cvec((float[]){2,0}, 1);
    complex32_t* B = cvec((float[]){3,4}, 1);
    complex32_t* C = cvec((float[]){1,1}, 1);
    cgemm('N', 'N', 1, 1, 1, alpha, A, 1, B, 1, beta, C, 1, 'n');
    assert_ulong(C[0].real.v, ==, 0x00000000u);                 // 0
    assert_ulong(C[0].imag.v, ==, 0x41800000u);                 // 16
    free(A); free(B); free(C);
    return MUNIT_OK;
}

//  B stored with a padding column (ldb=3 != ldc): exercises ldb indexing.
//  Real-only A=[[1,2],[3,4]], B=[[5,6],[7,8]] -> C=[[19,22],[43,50]].
MunitResult test_cgemm_ldb(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE;
    const complex32_t beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,0, 2,0, 3,0, 4,0}, 4);
    complex32_t* B = cvec((float[]){5,0, 6,0, 99,0,
                                    7,0, 8,0, 99,0}, 6);        // logical 2x2, ldb=3
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0}, 4);
    cgemm('N', 'N', 2, 2, 2, alpha, A, 2, B, 3, beta, C, 2, 'n');
    complex32_t* R = cvec((float[]){19,0, 22,0, 43,0, 50,0}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  Double-precision canonical vector (same matrices as test_cgemm_basic).
MunitResult test_zgemm_basic(const MunitParameter params[], void* u) {
    const complex64_t alpha = {{ 0x3ff0000000000000ull }, { 0 }};  // 1
    const complex64_t beta  = {{ 0 }, { 0 }};                       // 0
    complex64_t A[4] = {{{ 0x3ff0000000000000ull }, { 0x3ff0000000000000ull }},  // 1+1i
                        {{ 0x4000000000000000ull }, { 0 }},                       // 2
                        {{ 0 }, { 0 }},                                           // 0
                        {{ 0x3ff0000000000000ull }, { 0 }}};                      // 1
    complex64_t B[4] = {{{ 0x3ff0000000000000ull }, { 0 }},                       // 1
                        {{ 0 }, { 0 }},                                           // 0
                        {{ 0 }, { 0 }},                                           // 0
                        {{ 0x3ff0000000000000ull }, { 0x3ff0000000000000ull }}};  // 1+1i
    complex64_t C[4] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    zgemm('N', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    //  R = [[1+1i, 2+2i],[0, 1+1i]]
    assert_ullong(C[0].real.v, ==, 0x3ff0000000000000ull); assert_ullong(C[0].imag.v, ==, 0x3ff0000000000000ull);
    assert_ullong(C[1].real.v, ==, 0x4000000000000000ull); assert_ullong(C[1].imag.v, ==, 0x4000000000000000ull);
    assert_ullong(C[2].real.v, ==, 0x0ull);                assert_ullong(C[2].imag.v, ==, 0x0ull);
    assert_ullong(C[3].real.v, ==, 0x3ff0000000000000ull); assert_ullong(C[3].imag.v, ==, 0x3ff0000000000000ull);
    return MUNIT_OK;
}

//  Double-precision conjugate transpose (mirrors test_cgemm_conjtrans). A stored
//  [[1+1i, 2+1i],[3-1i, 0+2i]] -> op(A)=conj(A^T)=[[1-1i, 3+1i],[2-1i, 0-2i]].
MunitResult test_zgemm_conjtrans(const MunitParameter params[], void* u) {
    const complex64_t alpha = {{ 0x3ff0000000000000ull }, { 0 }};  // 1
    const complex64_t beta  = {{ 0 }, { 0 }};                       // 0
    const uint64_t P1 = 0x3ff0000000000000ull;   // 1
    const uint64_t N1 = 0xbff0000000000000ull;   // -1
    const uint64_t P2 = 0x4000000000000000ull;   // 2
    const uint64_t N2 = 0xc000000000000000ull;   // -2
    const uint64_t P3 = 0x4008000000000000ull;   // 3
    complex64_t A[4] = {{{ P1 }, { P1 }}, {{ P2 }, { P1 }},
                        {{ P3 }, { N1 }}, {{ 0 }, { P2 }}};        // [[1+1i,2+1i],[3-1i,0+2i]]
    complex64_t B[4] = {{{ P1 }, { 0 }}, {{ 0 }, { 0 }},
                        {{ 0 }, { 0 }}, {{ P1 }, { 0 }}};          // identity
    complex64_t C[4] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    zgemm('C', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    //  R = [[1-1i, 3+1i],[2-1i, 0-2i]]
    assert_ullong(C[0].real.v, ==, P1); assert_ullong(C[0].imag.v, ==, N1);
    assert_ullong(C[1].real.v, ==, P3); assert_ullong(C[1].imag.v, ==, P1);
    assert_ullong(C[2].real.v, ==, P2); assert_ullong(C[2].imag.v, ==, N1);
    assert_ullong(C[3].real.v, ==, 0x0ull); assert_ullong(C[3].imag.v, ==, N2);
    return MUNIT_OK;
}
