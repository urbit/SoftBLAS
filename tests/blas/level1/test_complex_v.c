#include "test.h"

//  Quad-precision complex (v) Level-1 routines. Exact small integers; the low
//  word is 0, so tests assert the high word v[1]. quad hi: 1=3fff.. 2=4000..
//  3=40008.. 4=40010.. 5=40014.. 6=40018.. 10=40024.. 11=40026.. -1=bfff..
//  -2=c000.. -3=c0008..
#define QH(hi) {{ 0x0ull, (hi) }}            // quad value from high word
#define ONE_HI   0x3fff000000000000ull
#define TWO_HI   0x4000000000000000ull
#define THREE_HI 0x4000800000000000ull
#define FOUR_HI  0x4001000000000000ull
#define FIVE_HI  0x4001400000000000ull
#define SIX_HI   0x4001800000000000ull
#define TEN_HI   0x4002400000000000ull
#define ELE_HI   0x4002600000000000ull
#define NEG2_HI  0xc000000000000000ull
#define NEG1_HI  0xbfff000000000000ull
#define NEG3_HI  0xc000800000000000ull

MunitResult test_vaxpy_basic(const MunitParameter params[], void* u) {
    const complex128_t CA = { QH(ONE_HI), QH(0x0ull) };               // 1
    complex128_t CX[1] = {{ QH(ONE_HI), QH(TWO_HI) }};               // 1 + 2i
    complex128_t CY[1] = {{ QH(THREE_HI), QH(FOUR_HI) }};            // 3 + 4i
    vaxpy(1, CA, CX, 1, CY, 1, 'n');                                 // 4 + 6i
    assert_ullong(CY[0].real.v[1], ==, FOUR_HI);
    assert_ullong(CY[0].imag.v[1], ==, SIX_HI);
    return MUNIT_OK;
}
MunitResult test_vcopy_basic(const MunitParameter params[], void* u) {
    complex128_t CX[1] = {{ QH(ONE_HI), QH(TWO_HI) }};
    complex128_t CY[1] = {{ QH(0x0ull), QH(0x0ull) }};
    vcopy(1, CX, 1, CY, 1, 'n');
    assert_ullong(CY[0].real.v[1], ==, ONE_HI);
    assert_ullong(CY[0].imag.v[1], ==, TWO_HI);
    return MUNIT_OK;
}
MunitResult test_vdotc_basic(const MunitParameter params[], void* u) {
    complex128_t CX[1] = {{ QH(ONE_HI), QH(TWO_HI) }};               // 1 + 2i
    complex128_t CY[1] = {{ QH(THREE_HI), QH(FOUR_HI) }};            // 3 + 4i
    complex128_t r = vdotc(1, CX, 1, CY, 1, 'n');                    // 11 - 2i
    assert_ullong(r.real.v[1], ==, ELE_HI);
    assert_ullong(r.imag.v[1], ==, NEG2_HI);
    return MUNIT_OK;
}
MunitResult test_vscal_basic(const MunitParameter params[], void* u) {
    const complex128_t CA = { QH(TWO_HI), QH(0x0ull) };              // 2
    complex128_t CX[1] = {{ QH(ONE_HI), QH(ONE_HI) }};              // 1 + 1i
    vscal(1, CA, CX, 1, 'n');                                        // 2 + 2i
    assert_ullong(CX[0].real.v[1], ==, TWO_HI);
    assert_ullong(CX[0].imag.v[1], ==, TWO_HI);
    return MUNIT_OK;
}
MunitResult test_vswap_basic(const MunitParameter params[], void* u) {
    complex128_t CX[1] = {{ QH(ONE_HI), QH(0x0ull) }};
    complex128_t CY[1] = {{ QH(TWO_HI), QH(0x0ull) }};
    vswap(1, CX, 1, CY, 1, 'n');
    assert_ullong(CX[0].real.v[1], ==, TWO_HI);
    assert_ullong(CY[0].real.v[1], ==, ONE_HI);
    return MUNIT_OK;
}
MunitResult test_vqrot_basic(const MunitParameter params[], void* u) {
    complex128_t CX[1] = {{ QH(ONE_HI), QH(0x0ull) }};
    complex128_t CY[1] = {{ QH(TWO_HI), QH(0x0ull) }};
    const complex128_t c = { QH(0x0ull), QH(0x0ull) }, s = { QH(ONE_HI), QH(0x0ull) };
    vqrot(1, CX, 1, CY, 1, c, s, 'n');
    assert_ullong(CX[0].real.v[1], ==, TWO_HI);    // 2
    assert_ullong(CY[0].real.v[1], ==, NEG1_HI);   // -1
    return MUNIT_OK;
}
MunitResult test_qvasum_basic(const MunitParameter params[], void* u) {
    complex128_t CX[2] = {{ QH(ONE_HI), QH(NEG2_HI) }, { QH(NEG3_HI), QH(FOUR_HI) }};  // 1,-2,-3,4
    assert_ullong(qvasum(2, CX, 1, 'n').v[1], ==, TEN_HI);                             // 10
    return MUNIT_OK;
}
MunitResult test_qvnrm2_basic(const MunitParameter params[], void* u) {
    complex128_t CX[1] = {{ QH(THREE_HI), QH(FOUR_HI) }};   // 3 + 4i
    assert_ullong(qvnrm2(1, CX, 1, 'n').v[1], ==, FIVE_HI); // 5
    return MUNIT_OK;
}
MunitResult test_ivamax_basic(const MunitParameter params[], void* u) {
    complex128_t CX[3] = {{ QH(ONE_HI), QH(0x0ull) }, { QH(0x0ull), QH(THREE_HI) }, { QH(TWO_HI), QH(0x0ull) }}; // 1,3,2
    assert_ulong(ivamax(3, CX, 1, 'n'), ==, 1u);
    return MUNIT_OK;
}
