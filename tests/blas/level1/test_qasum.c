#include "test.h"

MunitResult test_qasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = SB_REAL128U_ZERO, .lo = SB_REAL128L_ZERO},
            {.hi = SB_REAL128U_ZERO, .lo = SB_REAL128L_ZERO}},
        2);
    float128_t Q = (float128_t) qasum(2, (float128_t*)QX, 1);
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

    float128_t Q = (float128_t) qasum(5, (float128_t*)QX, 1);
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

    float128_t Q = (float128_t) qasum(5, (float128_t*)QX, 2);
    float128_t R = {0x0000000000000000, 0x4001400000000000};  // 5.0q

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}
