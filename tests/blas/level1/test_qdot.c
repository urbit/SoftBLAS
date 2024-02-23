#include "test.h"

MunitResult test_qdot_0(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        5);

    float128_t QZ = qdot(5, QX, 1, QY, 1);
    float128_t RZ = {0x0000000000000000, 0x0000000000000000};

    assert_ullong(QZ.v[0], ==, RZ.v[0]);
    assert_ullong(QZ.v[1], ==, RZ.v[1]);

    free(QX);
    free(QY);

    return MUNIT_OK;
}

MunitResult test_qdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000}},
        5);

    float128_t QZ = qdot(5, QX, 1, QY, 1);
    float128_t RZ = {0x0000000000000000, 0x4004180000000000};

    assert_ullong(QZ.v[0], ==, RZ.v[0]);
    assert_ullong(QZ.v[1], ==, RZ.v[1]);

    free(QX);
    free(QY);

    return MUNIT_OK;
}

MunitResult test_qdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
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

    float128_t QZ = qdot(5, QX, 1, QY, 2);
    float128_t RZ = {0x0000000000000000, 0x4002e00000000000};

    assert_ullong(QZ.v[0], ==, RZ.v[0]);
    assert_ullong(QZ.v[1], ==, RZ.v[1]);

    free(QX);
    free(QY);

    return MUNIT_OK;
}

MunitResult test_qdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000}},
        5);

    float128_t QZ = qdot(5, QX, 1, QY, -1);
    float128_t RZ = {0x0000000000000000, 0x4004b80000000000};

    assert_ullong(QZ.v[0], ==, RZ.v[0]);
    assert_ullong(QZ.v[1], ==, RZ.v[1]);

    free(QX);
    free(QY);

    return MUNIT_OK;
}