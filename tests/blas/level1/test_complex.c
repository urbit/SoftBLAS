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
