#include "test.h"

MunitResult test_qasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = SB_REAL128U_ZERO, .lo = SB_REAL128L_ZERO},
            {.hi = SB_REAL128U_ZERO, .lo = SB_REAL128L_ZERO}},
        2);
    float128_t Q = (float128_t) qasum(2, (float128_t*)QX, 1, 'n');
    float128_t R = {SB_REAL128L_ZERO, SB_REAL128U_ZERO};

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}

MunitResult test_qasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    float128_t Q = (float128_t) qasum(5, (float128_t*)QX, 1, 'n');
    float128_t R = {0x0000000000000000, 0x4002e00000000000};  // 15.0q

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}

MunitResult test_qasum_stride(const MunitParameter params[],
                              void *user_data) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000400000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000}},
        9);

    float128_t Q = (float128_t) qasum(5, (float128_t*)QX, 2, 'n');
    float128_t R = {0x0000000000000000, 0x4001400000000000};  // 5.0q

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}

//  Regression for the f128_abs bug. -pi has bit 63 set in its low word; the
//  old macro cleared that bit (it masked both 64-bit words instead of only the
//  sign word v[1]), silently corrupting the mantissa. |-pi| must equal +pi.
MunitResult test_qasum_negpi(const MunitParameter params[],
                             void* user_data_or_fixture) {
    // -pi in IEEE binary128: sign bit set in the high word; low word unchanged.
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0xc000921fb54442d1, .lo = 0x8469898cc51701b8}},
        1);

    float128_t Q = (float128_t) qasum(1, (float128_t*)QX, 1, 'n');

    // +pi: v[0] = low word (sign-bit untouched), v[1] = high word.
    float128_t R = {0x8469898cc51701b8, 0x4000921fb54442d1};

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);
    return MUNIT_OK;
}

//  Pins the float128_t storage layout every q* test depends on: the sign and
//  exponent live in the high word v[1]; v[0] is the low mantissa word.
MunitResult test_qasum_layout(const MunitParameter params[],
                              void* user_data_or_fixture) {
    //  1.0 in IEEE binary128, written {lo, hi}.
    float128_t one = {{ 0x0000000000000000, 0x3fff000000000000 }};
    assert_ullong(one.v[0], ==, 0x0000000000000000);  // low mantissa word
    assert_ullong(one.v[1], ==, 0x3fff000000000000);  // sign+exponent word

    //  A routine that reads the pair agrees: |1.0| == 1.0, bit-for-bit.
    float128_t r = qasum(1, &one, 1, 'n');
    assert_ullong(r.v[1], ==, 0x3fff000000000000);
    assert_ullong(r.v[0], ==, 0x0000000000000000);
    return MUNIT_OK;
}
