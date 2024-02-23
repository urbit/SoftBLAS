#include "test.h"

MunitResult test_iqamax_0(const MunitParameter params[],
                          void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        3);

    uint64_t I = iqamax(3, QX, 1);
    uint64_t R = 0;

    assert_int(I, ==, R);
    
    free(QX);

    return MUNIT_OK;
}

MunitResult test_iqamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    uint64_t I = iqamax(5, QX, 1);
    uint64_t R = 4;

    assert_int(I, ==, R);
    
    free(QX);

    return MUNIT_OK;
}

MunitResult test_iqamax_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    uint64_t I = iqamax(2, QX, 3);
    // note that iqamax counts against the stride index, not the overall index
    uint64_t R = 1;

    assert_int(I, ==, R);
    
    free(QX);

    return MUNIT_OK;
}

MunitResult test_iqamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000}},
        5);

    uint64_t I = iqamax(5, QX, 1);
    uint64_t R = 2;

    assert_int(I, ==, R);
    
    free(QX);

    return MUNIT_OK;
}
