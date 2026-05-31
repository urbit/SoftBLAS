#include "test.h"

//  Half-precision complex (i) Level-1 routines. Small integer values keep every
//  result exact. half: 1=3c00 2=4000 3=4200 4=4400 5=4500 6=4600 10=4900
//  11=4980 -1=bc00 -2=c000.

MunitResult test_iaxpy_basic(const MunitParameter params[], void* u) {
    const complex16_t CA = {{ 0x3c00 }, { 0x0000 }};            // 1 + 0i
    complex16_t CX[1] = {{{ 0x3c00 }, { 0x4000 }}};            // 1 + 2i
    complex16_t CY[1] = {{{ 0x4200 }, { 0x4400 }}};            // 3 + 4i
    iaxpy(1, CA, CX, 1, CY, 1, 'n');                          // -> 4 + 6i
    assert_ushort(CY[0].real.v, ==, 0x4400);                 // 4
    assert_ushort(CY[0].imag.v, ==, 0x4600);                 // 6
    return MUNIT_OK;
}
MunitResult test_icopy_basic(const MunitParameter params[], void* u) {
    complex16_t CX[1] = {{{ 0x3c00 }, { 0x4000 }}};           // 1 + 2i
    complex16_t CY[1] = {{{ 0 }, { 0 }}};
    icopy(1, CX, 1, CY, 1, 'n');
    assert_ushort(CY[0].real.v, ==, 0x3c00);
    assert_ushort(CY[0].imag.v, ==, 0x4000);
    return MUNIT_OK;
}
MunitResult test_idotc_basic(const MunitParameter params[], void* u) {
    complex16_t CX[1] = {{{ 0x3c00 }, { 0x4000 }}};           // 1 + 2i
    complex16_t CY[1] = {{{ 0x4200 }, { 0x4400 }}};           // 3 + 4i
    complex16_t r = idotc(1, CX, 1, CY, 1, 'n');             // conj(x)*y = 11 - 2i
    assert_ushort(r.real.v, ==, 0x4980);                     // 11
    assert_ushort(r.imag.v, ==, 0xc000);                     // -2
    return MUNIT_OK;
}
MunitResult test_iscal_basic(const MunitParameter params[], void* u) {
    const complex16_t CA = {{ 0x4000 }, { 0x0000 }};          // 2 + 0i
    complex16_t CX[1] = {{{ 0x3c00 }, { 0x3c00 }}};          // 1 + 1i
    iscal(1, CA, CX, 1, 'n');                                // 2 + 2i
    assert_ushort(CX[0].real.v, ==, 0x4000);
    assert_ushort(CX[0].imag.v, ==, 0x4000);
    return MUNIT_OK;
}
MunitResult test_iswap_basic(const MunitParameter params[], void* u) {
    complex16_t CX[1] = {{{ 0x3c00 }, { 0 }}};               // 1
    complex16_t CY[1] = {{{ 0x4000 }, { 0 }}};               // 2
    iswap(1, CX, 1, CY, 1, 'n');
    assert_ushort(CX[0].real.v, ==, 0x4000);                 // 2
    assert_ushort(CY[0].real.v, ==, 0x3c00);                 // 1
    return MUNIT_OK;
}
MunitResult test_ihrot_basic(const MunitParameter params[], void* u) {
    complex16_t CX[1] = {{{ 0x3c00 }, { 0 }}};               // 1
    complex16_t CY[1] = {{{ 0x4000 }, { 0 }}};               // 2
    const complex16_t c = {{ 0 }, { 0 }}, s = {{ 0x3c00 }, { 0 }};  // 0, 1
    ihrot(1, CX, 1, CY, 1, c, s, 'n');
    assert_ushort(CX[0].real.v, ==, 0x4000);                 // 2
    assert_ushort(CY[0].real.v, ==, 0xbc00);                 // -1
    return MUNIT_OK;
}
MunitResult test_hiasum_basic(const MunitParameter params[], void* u) {
    complex16_t CX[2] = {{{ 0x3c00 }, { 0xc000 }}, {{ 0xc200 }, { 0x4400 }}};  // 1,-2,-3,4
    assert_ushort(hiasum(2, CX, 1, 'n').v, ==, 0x4900);      // |1|+|2|+|3|+|4| = 10
    return MUNIT_OK;
}
MunitResult test_hinrm2_basic(const MunitParameter params[], void* u) {
    complex16_t CX[1] = {{{ 0x4200 }, { 0x4400 }}};          // 3 + 4i
    assert_ushort(hinrm2(1, CX, 1, 'n').v, ==, 0x4500);      // 5
    return MUNIT_OK;
}
MunitResult test_iiamax_basic(const MunitParameter params[], void* u) {
    complex16_t CX[3] = {{{ 0x3c00 }, { 0 }}, {{ 0 }, { 0x4200 }}, {{ 0x4000 }, { 0 }}};  // mags 1,3,2
    assert_ulong(iiamax(3, CX, 1, 'n'), ==, 1u);
    return MUNIT_OK;
}
//  NaN in either component canonicalizes both to HALFNAN (nan_unify_i).
MunitResult test_iaxpy_nan(const MunitParameter params[], void* u) {
    const complex16_t CA = {{ 0x3c00 }, { 0 }};               // 1
    complex16_t CX[1] = {{{ 0x7e01 }, { 0 }}};               // NaN(non-canon) + 0i
    complex16_t CY[1] = {{{ 0 }, { 0 }}};
    iaxpy(1, CA, CX, 1, CY, 1, 'n');
    assert_ushort(CY[0].real.v, ==, (uint16_t)HALFNAN);
    assert_ushort(CY[0].imag.v, ==, (uint16_t)HALFNAN);
    return MUNIT_OK;
}
//  Negative stride reverses the indexed vector (offset start (-N+1)*incX).
MunitResult test_iaxpy_negstride(const MunitParameter params[], void* u) {
    const complex16_t CA = {{ 0x3c00 }, { 0 }};               // 1
    complex16_t CX[3] = {{{ 0x3c00 }, { 0 }}, {{ 0x4000 }, { 0 }}, {{ 0x4200 }, { 0 }}};  // 1,2,3
    complex16_t CY[3] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    iaxpy(3, CA, CX, -1, CY, 1, 'n');                         // CY += reverse(CX)
    assert_ushort(CY[0].real.v, ==, 0x4200);                 // 3
    assert_ushort(CY[2].real.v, ==, 0x3c00);                 // 1
    return MUNIT_OK;
}
MunitResult test_icopy_negstride(const MunitParameter params[], void* u) {
    complex16_t CX[3] = {{{ 0x3c00 }, { 0 }}, {{ 0x4000 }, { 0 }}, {{ 0x4200 }, { 0 }}};  // 1,2,3
    complex16_t CY[3] = {{{ 0 }, { 0 }}, {{ 0 }, { 0 }}, {{ 0 }, { 0 }}};
    icopy(3, CX, -1, CY, 1, 'n');                            // CY = reverse(CX)
    assert_ushort(CY[0].real.v, ==, 0x4200);                 // 3
    assert_ushort(CY[2].real.v, ==, 0x3c00);                 // 1
    return MUNIT_OK;
}
