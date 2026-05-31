#include "test.h"

//  Complex (single-precision) Level-1 routines.

MunitResult test_caxpy_basic(const MunitParameter params[], void* u) {
    const complex32_t CA = {{ 0x3f800000 }, { 0x00000000 }};   // 1 + 0i
    complex32_t* CX = cvec((float[]){1.0f, 2.0f}, 1);          // 1 + 2i
    complex32_t* CY = cvec((float[]){3.0f, 4.0f}, 1);          // 3 + 4i
    caxpy(1, CA, CX, 1, CY, 1, 'n');
    assert_ulong(CY[0].real.v, ==, 0x40800000u);              // 4
    assert_ulong(CY[0].imag.v, ==, 0x40c00000u);              // 6i
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_ccopy_basic(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1.0f, 2.0f}, 1);
    complex32_t* CY = cvec((float[]){0.0f, 0.0f}, 1);
    ccopy(1, CX, 1, CY, 1, 'n');
    assert_ulong(CY[0].real.v, ==, 0x3f800000u);              // 1
    assert_ulong(CY[0].imag.v, ==, 0x40000000u);              // 2i
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_cdotc_conj(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1.0f, 2.0f}, 1);          // 1 + 2i
    complex32_t* CY = cvec((float[]){3.0f, 4.0f}, 1);          // 3 + 4i
    complex32_t r = cdotc(1, CX, 1, CY, 1, 'n');              // conj(x)*y = 11 - 2i
    assert_ulong(r.real.v, ==, 0x41300000u);                  // 11
    assert_ulong(r.imag.v, ==, 0xc0000000u);                  // -2i
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_scasum_basic(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1.0f, -2.0f, -3.0f, 4.0f}, 2);  // |1|+|2|+|3|+|4|
    float32_t r = scasum(2, CX, 1, 'n');
    assert_ulong(r.v, ==, 0x41200000u);                      // 10
    free(CX);
    return MUNIT_OK;
}
MunitResult test_scnrm2_basic(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){3.0f, 4.0f}, 1);         // |3 + 4i| = 5
    float32_t r = scnrm2(1, CX, 1, 'n');
    assert_ulong(r.v, ==, 0x40a00000u);                      // 5
    free(CX);
    return MUNIT_OK;
}
MunitResult test_csrot_basic(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1.0f, 0.0f}, 1);
    complex32_t* CY = cvec((float[]){2.0f, 0.0f}, 1);
    const complex32_t c = {{ 0x00000000 }, { 0x00000000 }};  // 0
    const complex32_t s = {{ 0x3f800000 }, { 0x00000000 }};  // 1
    csrot(1, CX, 1, CY, 1, c, s, 'n');                       // 90-degree-like swap
    assert_ulong(CX[0].real.v, ==, 0x40000000u);             // 2
    assert_ulong(CY[0].real.v, ==, 0xbf800000u);             // -1
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_cscal_basic(const MunitParameter params[], void* u) {
    const complex32_t CA = {{ 0x40000000 }, { 0x00000000 }};  // 2 + 0i
    complex32_t* CX = cvec((float[]){1.0f, 1.0f}, 1);         // 1 + 1i
    cscal(1, CA, CX, 1, 'n');                                // (2)(1+1i) = 2 + 2i
    assert_ulong(CX[0].real.v, ==, 0x40000000u);             // 2
    assert_ulong(CX[0].imag.v, ==, 0x40000000u);             // 2i
    free(CX);
    return MUNIT_OK;
}
MunitResult test_cswap_basic(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1.0f, 0.0f}, 1);
    complex32_t* CY = cvec((float[]){2.0f, 0.0f}, 1);
    cswap(1, CX, 1, CY, 1, 'n');
    assert_ulong(CX[0].real.v, ==, 0x40000000u);             // 2
    assert_ulong(CY[0].real.v, ==, 0x3f800000u);             // 1
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_icamax_basic(const MunitParameter params[], void* u) {
    // |Re|+|Im| magnitudes: 1, 3, 2 -> max at 0-based index 1
    complex32_t* CX = cvec((float[]){1.0f, 0.0f, 0.0f, 3.0f, 2.0f, 0.0f}, 3);
    assert_ulong(icamax(3, CX, 1, 'n'), ==, 1u);
    assert_ulong(icamax(1, CX, 1, 'n'), ==, 0u);             // N==1 -> 0
    free(CX);
    return MUNIT_OK;
}

//  C8: complex NaN canonicalization. nan_unify_c canonicalizes BOTH components
//  to SINGNAN when either is a (non-canonical) NaN, so the determinism
//  guarantee holds for the complex routines too.
#define CNCNAN 0x7fc00001u   // non-canonical single NaN
MunitResult test_caxpy_nan(const MunitParameter params[], void* u) {
    const complex32_t CA = {{ 0x3f800000 }, { 0x00000000 }};   // 1 + 0i
    complex32_t CX[1] = {{{ CNCNAN }, { 0x00000000 }}};        // NaN + 0i
    complex32_t CY[1] = {{{ 0x00000000 }, { 0x00000000 }}};
    caxpy(1, CA, CX, 1, CY, 1, 'n');
    assert_ulong(CY[0].real.v, ==, (uint32_t)SINGNAN);
    assert_ulong(CY[0].imag.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_cscal_nan(const MunitParameter params[], void* u) {
    const complex32_t CA = {{ 0x40000000 }, { 0x00000000 }};   // 2 + 0i
    complex32_t CX[1] = {{{ CNCNAN }, { 0x3f800000 }}};        // NaN + 1i
    cscal(1, CA, CX, 1, 'n');
    assert_ulong(CX[0].real.v, ==, (uint32_t)SINGNAN);
    assert_ulong(CX[0].imag.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_ccopy_nan(const MunitParameter params[], void* u) {
    complex32_t CX[1] = {{{ 0x00000000 }, { CNCNAN }}};        // 0 + NaN i
    complex32_t CY[1] = {{{ 0x00000000 }, { 0x00000000 }}};
    ccopy(1, CX, 1, CY, 1, 'n');
    assert_ulong(CY[0].real.v, ==, (uint32_t)SINGNAN);        // copied NaN canonicalized
    assert_ulong(CY[0].imag.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}

//  C8: negative-stride coverage for the complex family. incX=-1 reverses the
//  indexed vector (offset start (-N+1)*incX).
MunitResult test_caxpy_negstride(const MunitParameter params[], void* u) {
    const complex32_t CA = {{ 0x3f800000 }, { 0x00000000 }};   // 1 + 0i
    complex32_t* CX = cvec((float[]){1,0, 2,0, 3,0}, 3);
    complex32_t* CY = cvec((float[]){0,0, 0,0, 0,0}, 3);
    caxpy(3, CA, CX, -1, CY, 1, 'n');                          // CY += reverse(CX)
    assert_ulong(CY[0].real.v, ==, 0x40400000u);              // 3
    assert_ulong(CY[1].real.v, ==, 0x40000000u);              // 2
    assert_ulong(CY[2].real.v, ==, 0x3f800000u);              // 1
    free(CX); free(CY);
    return MUNIT_OK;
}
MunitResult test_ccopy_negstride(const MunitParameter params[], void* u) {
    complex32_t* CX = cvec((float[]){1,0, 2,0, 3,0}, 3);
    complex32_t* CY = cvec((float[]){0,0, 0,0, 0,0}, 3);
    ccopy(3, CX, -1, CY, 1, 'n');                             // CY = reverse(CX)
    assert_ulong(CY[0].real.v, ==, 0x40400000u);              // 3
    assert_ulong(CY[2].real.v, ==, 0x3f800000u);              // 1
    free(CX); free(CY);
    return MUNIT_OK;
}

//  B5: the complex asum/nrm2 reductions (all four precisions) must reject a
//  stride where (N-1)*incX overflows uint64_t, returning canonical 0 instead of
//  reading off the (single-element) buffer. ASan-clean.
MunitResult test_complex_huge_stride(const MunitParameter params[], void* u) {
    const uint64_t H = 0x1000000000ull;  // 2^36; (H-1)*H wraps uint64_t
    complex32_t  c1[1] = {{{ 0x3f800000u }, { 0 }}};
    complex16_t  i1[1] = {{{ 0x3c00 }, { 0 }}};
    complex64_t  z1[1] = {{{ 0x3ff0000000000000ull }, { 0 }}};
    complex128_t v1[1] = {{{{ 0, 0x3fff000000000000ull }}, {{ 0, 0 }}}};
    assert_ulong(scasum(H, c1, H, 'n').v, ==, 0x0u);
    assert_ulong(scnrm2(H, c1, H, 'n').v, ==, 0x0u);
    assert_ushort(hiasum(H, i1, H, 'n').v, ==, 0x0);
    assert_ushort(hinrm2(H, i1, H, 'n').v, ==, 0x0);
    assert_ullong(dzasum(H, z1, H, 'n').v, ==, 0x0ull);
    assert_ullong(dznrm2(H, z1, H, 'n').v, ==, 0x0ull);
    assert_ullong(qvasum(H, v1, H, 'n').v[1], ==, 0x0ull);
    assert_ullong(qvnrm2(H, v1, H, 'n').v[1], ==, 0x0ull);
    return MUNIT_OK;
}
