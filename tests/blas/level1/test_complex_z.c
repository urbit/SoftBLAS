#include "test.h"

//  Double-precision complex (z) Level-1 routines. Exact small integers.
//  double: 1=3ff0.. 2=4000.. 3=4008.. 4=4010.. 5=4014.. 6=4018.. 10=4024..
//  11=4026.. -1=bff0.. -2=c000.. -3=c008..

MunitResult test_zaxpy_basic(const MunitParameter params[], void* u) {
    const complex64_t CA = {{ 0x3ff0000000000000ull }, { 0 }};                    // 1
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0x4000000000000000ull }}}; // 1 + 2i
    complex64_t CY[1] = {{{ 0x4008000000000000ull }, { 0x4010000000000000ull }}}; // 3 + 4i
    zaxpy(1, CA, CX, 1, CY, 1, 'n');                                              // 4 + 6i
    assert_ullong(CY[0].real.v, ==, 0x4010000000000000ull);                      // 4
    assert_ullong(CY[0].imag.v, ==, 0x4018000000000000ull);                      // 6
    return MUNIT_OK;
}
MunitResult test_zcopy_basic(const MunitParameter params[], void* u) {
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0x4000000000000000ull }}};
    complex64_t CY[1] = {{{ 0 }, { 0 }}};
    zcopy(1, CX, 1, CY, 1, 'n');
    assert_ullong(CY[0].real.v, ==, 0x3ff0000000000000ull);
    assert_ullong(CY[0].imag.v, ==, 0x4000000000000000ull);
    return MUNIT_OK;
}
MunitResult test_zdotc_basic(const MunitParameter params[], void* u) {
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0x4000000000000000ull }}}; // 1 + 2i
    complex64_t CY[1] = {{{ 0x4008000000000000ull }, { 0x4010000000000000ull }}}; // 3 + 4i
    complex64_t r = zdotc(1, CX, 1, CY, 1, 'n');                                  // 11 - 2i
    assert_ullong(r.real.v, ==, 0x4026000000000000ull);                          // 11
    assert_ullong(r.imag.v, ==, 0xc000000000000000ull);                          // -2
    return MUNIT_OK;
}
MunitResult test_zscal_basic(const MunitParameter params[], void* u) {
    const complex64_t CA = {{ 0x4000000000000000ull }, { 0 }};                    // 2
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0x3ff0000000000000ull }}}; // 1 + 1i
    zscal(1, CA, CX, 1, 'n');                                                     // 2 + 2i
    assert_ullong(CX[0].real.v, ==, 0x4000000000000000ull);
    assert_ullong(CX[0].imag.v, ==, 0x4000000000000000ull);
    return MUNIT_OK;
}
MunitResult test_zswap_basic(const MunitParameter params[], void* u) {
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0 }}};
    complex64_t CY[1] = {{{ 0x4000000000000000ull }, { 0 }}};
    zswap(1, CX, 1, CY, 1, 'n');
    assert_ullong(CX[0].real.v, ==, 0x4000000000000000ull);
    assert_ullong(CY[0].real.v, ==, 0x3ff0000000000000ull);
    return MUNIT_OK;
}
MunitResult test_zdrot_basic(const MunitParameter params[], void* u) {
    complex64_t CX[1] = {{{ 0x3ff0000000000000ull }, { 0 }}};
    complex64_t CY[1] = {{{ 0x4000000000000000ull }, { 0 }}};
    const complex64_t c = {{ 0 }, { 0 }}, s = {{ 0x3ff0000000000000ull }, { 0 }};
    zdrot(1, CX, 1, CY, 1, c, s, 'n');
    assert_ullong(CX[0].real.v, ==, 0x4000000000000000ull);   // 2
    assert_ullong(CY[0].real.v, ==, 0xbff0000000000000ull);   // -1
    return MUNIT_OK;
}
MunitResult test_dzasum_basic(const MunitParameter params[], void* u) {
    complex64_t CX[2] = {{{ 0x3ff0000000000000ull }, { 0xc000000000000000ull }},
                         {{ 0xc008000000000000ull }, { 0x4010000000000000ull }}}; // 1,-2,-3,4
    assert_ullong(dzasum(2, CX, 1, 'n').v, ==, 0x4024000000000000ull);            // 10
    return MUNIT_OK;
}
MunitResult test_dznrm2_basic(const MunitParameter params[], void* u) {
    complex64_t CX[1] = {{{ 0x4008000000000000ull }, { 0x4010000000000000ull }}}; // 3 + 4i
    assert_ullong(dznrm2(1, CX, 1, 'n').v, ==, 0x4014000000000000ull);            // 5
    return MUNIT_OK;
}
MunitResult test_izamax_basic(const MunitParameter params[], void* u) {
    complex64_t CX[3] = {{{ 0x3ff0000000000000ull }, { 0 }},
                         {{ 0 }, { 0x4008000000000000ull }},
                         {{ 0x4000000000000000ull }, { 0 }}};                     // mags 1,3,2
    assert_ulong(izamax(3, CX, 1, 'n'), ==, 1u);
    return MUNIT_OK;
}
