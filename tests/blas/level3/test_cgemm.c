#include "test.h"

//  Complex GEMM (half 'i' / single 'c' / double 'z' / quad 'v').
//  C <- alpha*op(A)*op(B) + beta*C, op = N/T/C. Row-major, leading dims
//  lda/ldb/ldc, mirroring sgemm. All test values are small integers so the
//  float results are exact and bit-checkable.

//  Quad helpers (high word holds the value; low word stays 0).
#define QH(hi) {{ 0x0ull, (hi) }}
#define Q_ONE 0x3fff000000000000ull
#define Q_TWO 0x4000000000000000ull

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

//  Half-precision canonical vector (same matrices as test_cgemm_basic).
MunitResult test_igemm_basic(const MunitParameter params[], void* u) {
    const complex16_t alpha = SB_COMPLEX16_ONE;
    const complex16_t beta = SB_COMPLEX16_ZERO;
    complex16_t A[4] = {{{ 0x3c00 }, { 0x3c00 }},  // 1+1i
                        {{ 0x4000 }, { 0 }},        // 2
                        {{ 0 }, { 0 }},             // 0
                        {{ 0x3c00 }, { 0 }}};       // 1
    complex16_t B[4] = {{{ 0x3c00 }, { 0 }},        // 1
                        {{ 0 }, { 0 }},             // 0
                        {{ 0 }, { 0 }},             // 0
                        {{ 0x3c00 }, { 0x3c00 }}};  // 1+1i
    complex16_t C[4] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    igemm('N', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    //  R = [[1+1i, 2+2i],[0, 1+1i]]
    assert_ushort(C[0].real.v, ==, 0x3c00); assert_ushort(C[0].imag.v, ==, 0x3c00);
    assert_ushort(C[1].real.v, ==, 0x4000); assert_ushort(C[1].imag.v, ==, 0x4000);
    assert_ushort(C[2].real.v, ==, 0x0);    assert_ushort(C[2].imag.v, ==, 0x0);
    assert_ushort(C[3].real.v, ==, 0x3c00); assert_ushort(C[3].imag.v, ==, 0x3c00);
    return MUNIT_OK;
}

//  Quad-precision canonical vector (same matrices as test_cgemm_basic).
MunitResult test_vgemm_basic(const MunitParameter params[], void* u) {
    const complex128_t alpha = { QH(Q_ONE), QH(0x0ull) };  // 1
    const complex128_t beta  = { QH(0x0ull), QH(0x0ull) };  // 0
    complex128_t A[4] = {{ QH(Q_ONE), QH(Q_ONE) },          // 1+1i
                         { QH(Q_TWO), QH(0x0ull) },         // 2
                         { QH(0x0ull), QH(0x0ull) },        // 0
                         { QH(Q_ONE), QH(0x0ull) }};        // 1
    complex128_t B[4] = {{ QH(Q_ONE), QH(0x0ull) },         // 1
                         { QH(0x0ull), QH(0x0ull) },        // 0
                         { QH(0x0ull), QH(0x0ull) },        // 0
                         { QH(Q_ONE), QH(Q_ONE) }};         // 1+1i
    complex128_t C[4] = {{ QH(0x0ull), QH(0x0ull) }, { QH(0x0ull), QH(0x0ull) },
                         { QH(0x0ull), QH(0x0ull) }, { QH(0x0ull), QH(0x0ull) }};
    vgemm('N', 'N', 2, 2, 2, alpha, A, 2, B, 2, beta, C, 2, 'n');
    //  R = [[1+1i, 2+2i],[0, 1+1i]]
    assert_ullong(C[0].real.v[1], ==, Q_ONE); assert_ullong(C[0].imag.v[1], ==, Q_ONE);
    assert_ullong(C[1].real.v[1], ==, Q_TWO); assert_ullong(C[1].imag.v[1], ==, Q_TWO);
    assert_ullong(C[2].real.v[1], ==, 0x0ull); assert_ullong(C[2].imag.v[1], ==, 0x0ull);
    assert_ullong(C[3].real.v[1], ==, Q_ONE); assert_ullong(C[3].imag.v[1], ==, Q_ONE);
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

//  3x3 single (N,N). A,B,C row-major, lda=ldb=ldc=3. Larger than the 2x2 cases
//  so a transposed-index bug in the triple loop would surface.
MunitResult test_cgemm_3x3(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,0, 0,1, 0,2, 1,1, 3,0, 1,0, 0,1, 2,2}, 9);
    complex32_t* B = cvec((float[]){1,0, 0,1, 2,0, 1,1, 2,0, 0,1, 0,1, 1,0, 1,1}, 9);
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0}, 9);
    cgemm('N', 'N', 3, 3, 3, alpha, A, 3, B, 3, beta, C, 3, 'n');
    complex32_t* R = cvec((float[]){2,3, 3,2, 1,5, 0,7, 3,2, 2,8, -2,3, 2,5, 1,4}, 9);
    for (uint64_t i = 0; i < 9; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  4x4 single conjugate transpose: C = conj(A^T) * B, lda=ldb=ldc=4. Exercises
//  the 'C' op and ldc row-stepping at a size where indexing mistakes show.
MunitResult test_cgemm_4x4_conjtrans(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,-1, 0,1, 1,0, 0,2, 1,1, 3,0, 2,-2,
                                    1,0, 0,1, 2,2, 1,1, 3,-1, 1,0, 0,-1, 2,0}, 16);
    complex32_t* B = cvec((float[]){1,0, 0,1, 2,0, 1,-1, 1,1, 2,0, 0,1, 0,1,
                                    0,1, 1,0, 1,1, 2,0, 1,0, 1,1, 0,1, 1,0}, 16);
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
                                    0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0}, 16);
    cgemm('C', 'N', 4, 4, 4, alpha, A, 4, B, 4, beta, C, 4, 'n');
    complex32_t* R = cvec((float[]){6,-1, 4,1, 4,2, 7,-1, 6,1, 2,0, 6,3, 5,-2,
                                    5,5, 8,-1, 3,1, 3,-1, 4,5, 7,6, 2,4, 3,-1}, 16);
    for (uint64_t i = 0; i < 16; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  3x3 double (N,N), same matrices as test_cgemm_3x3.
MunitResult test_zgemm_3x3(const MunitParameter params[], void* u) {
    const complex64_t alpha = {{ 0x3ff0000000000000ull }, { 0 }}, beta = {{ 0 }, { 0 }};
    complex64_t* A = zvec((double[]){1,1, 2,0, 0,1, 0,2, 1,1, 3,0, 1,0, 0,1, 2,2}, 9);
    complex64_t* B = zvec((double[]){1,0, 0,1, 2,0, 1,1, 2,0, 0,1, 0,1, 1,0, 1,1}, 9);
    complex64_t* C = zvec((double[]){0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0}, 9);
    zgemm('N', 'N', 3, 3, 3, alpha, A, 3, B, 3, beta, C, 3, 'n');
    complex64_t* R = zvec((double[]){2,3, 3,2, 1,5, 0,7, 3,2, 2,8, -2,3, 2,5, 1,4}, 9);
    for (uint64_t i = 0; i < 9; i++) {
        assert_ullong(C[i].real.v, ==, R[i].real.v);
        assert_ullong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(A); free(B); free(C); free(R);
    return MUNIT_OK;
}

//  ------------------------------------------------------------------------
//  Closed-form / torture cases. Still bit-exact: every value below is a
//  Gaussian integer (entries in {0, ±1, ±i, ...}) so no rounding occurs.
//  ------------------------------------------------------------------------

//  DFT-4 unitarity. The 4-point DFT matrix has entries in {1, i, -1, -i}
//  (4th roots of unity, all exactly representable), and F·Fᴴ = 4·I exactly.
//  Computed as cgemm('N','C', F, F): op(B)=conj(Fᵀ)=Fᴴ. A scrambled index or a
//  missing conjugate in the 'C' path breaks the clean 4·I.
MunitResult test_cgemm_dft4_unitary(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* F = cvec((float[]){1,0,  1,0,  1,0,  1,0,
                                    1,0,  0,1, -1,0,  0,-1,
                                    1,0, -1,0,  1,0, -1,0,
                                    1,0,  0,-1, -1,0,  0,1}, 16);
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
                                    0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0}, 16);
    cgemm('N', 'C', 4, 4, 4, alpha, F, 4, F, 4, beta, C, 4, 'n');
    complex32_t* R = cvec((float[]){4,0, 0,0, 0,0, 0,0,  0,0, 4,0, 0,0, 0,0,
                                    0,0, 0,0, 4,0, 0,0,  0,0, 0,0, 0,0, 4,0}, 16);
    for (uint64_t i = 0; i < 16; i++) {
        assert_ulong(C[i].real.v, ==, R[i].real.v);
        assert_ulong(C[i].imag.v, ==, R[i].imag.v);
    }
    free(F); free(C); free(R);
    return MUNIT_OK;
}

//  Pauli algebra: σx·σy = i·σz and σy·σy = I. Genuine complex products
//  (σy = [[0,-i],[i,0]]) with exact ±1/±i entries.
MunitResult test_cgemm_pauli(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* SX = cvec((float[]){0,0, 1,0, 1,0, 0,0}, 4);   // σx
    complex32_t* SY = cvec((float[]){0,0, 0,-1, 0,1, 0,0}, 4);  // σy
    //  σx·σy = i·σz = [[i,0],[0,-i]]
    complex32_t* C = cvec((float[]){0,0, 0,0, 0,0, 0,0}, 4);
    cgemm('N', 'N', 2, 2, 2, alpha, SX, 2, SY, 2, beta, C, 2, 'n');
    complex32_t* R1 = cvec((float[]){0,1, 0,0, 0,0, 0,-1}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C[i].real.v, ==, R1[i].real.v);
        assert_ulong(C[i].imag.v, ==, R1[i].imag.v);
    }
    //  σy·σy = I
    complex32_t* C2 = cvec((float[]){0,0, 0,0, 0,0, 0,0}, 4);
    cgemm('N', 'N', 2, 2, 2, alpha, SY, 2, SY, 2, beta, C2, 2, 'n');
    complex32_t* R2 = cvec((float[]){1,0, 0,0, 0,0, 1,0}, 4);
    for (uint64_t i = 0; i < 4; i++) {
        assert_ulong(C2[i].real.v, ==, R2[i].real.v);
        assert_ulong(C2[i].imag.v, ==, R2[i].imag.v);
    }
    free(SX); free(SY); free(C); free(R1); free(C2); free(R2);
    return MUNIT_OK;
}

//  Gram matrix invariant: G = Aᴴ·A is Hermitian (G[i][j] = conj(G[j][i]), real
//  diagonal) for ANY A. Checked structurally via cgemm('C','N',A,A) — exercises
//  the 'C' op without a hand-computed expected matrix. Integer A keeps G[i][j]
//  and conj(G[j][i]) (different summation orders) bit-equal; f32_eq folds ±0.
MunitResult test_cgemm_hermitian_gram(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t* A = cvec((float[]){1,1, 2,-1, 0,1, 1,0, 0,2, 1,1, 3,0, 2,-2,
                                    1,0, 0,1, 2,2, 1,1, 3,-1, 1,0, 0,-1, 2,0}, 16);
    complex32_t* G = cvec((float[]){0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
                                    0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0}, 16);
    cgemm('C', 'N', 4, 4, 4, alpha, A, 4, A, 4, beta, G, 4, 'n');
    for (uint64_t i = 0; i < 4; i++) {
        for (uint64_t j = 0; j < 4; j++) {
            complex32_t gij = G[i*4 + j];
            complex32_t gji = G[j*4 + i];
            assert_true(c32_eq(gij, c32_conj(gji)));  // i==j -> imag must be ±0
        }
    }
    free(A); free(G);
    return MUNIT_OK;
}

//  NaN canonicalization: a NaN anywhere in A taints C[0], and nan_unify_c
//  forces BOTH components to SINGNAN (mirrors test_sgemm_nan).
MunitResult test_cgemm_nan(const MunitParameter params[], void* u) {
    const complex32_t alpha = SB_COMPLEX32_ONE, beta = SB_COMPLEX32_ZERO;
    complex32_t A[1] = {{{ 0x7fc00001u }, { 0 }}};   // NaN(non-canon) + 0i
    complex32_t B[1] = {{{ 0x3f800000u }, { 0 }}};   // 1 + 0i
    complex32_t C[1] = {{{ 0 }, { 0 }}};
    cgemm('N', 'N', 1, 1, 1, alpha, A, 1, B, 1, beta, C, 1, 'n');
    assert_ulong(C[0].real.v, ==, (uint32_t)SINGNAN);
    assert_ulong(C[0].imag.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
