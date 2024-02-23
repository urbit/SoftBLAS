#include "test.h"

MunitResult test_qnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        2);

    float128_t Q = (float128_t) qnrm2(2, (float128_t*)QX, 1);
    float128_t R = {0x0000000000000000, 0x0000000000000000};

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);
    
    free(QX);

    return MUNIT_OK;
}

MunitResult test_qnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    float128_t Q = (float128_t) qnrm2(5, (float128_t*)QX, 1);
    float128_t R = {0x7f600f8fe05dd35e, 0x4001daa2fefaae1d};  // sqrt(55.0q)

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}

MunitResult test_qnrm2_stride(const MunitParameter params[],
                              void *user_data) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000}},
        9);

    float128_t Q = (float128_t) qnrm2(5, (float128_t*)QX, 2);
    float128_t R = {0x7c15f39cc0605cee, 0x40001e3779b97f4a};  // sqrt(5.0q)

    assert_ullong(Q.v[0], ==, R.v[0]);
    assert_ullong(Q.v[1], ==, R.v[1]);

    free(QX);

    return MUNIT_OK;
}
