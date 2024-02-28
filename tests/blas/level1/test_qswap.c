#include "test.h"

MunitResult test_qswap_two(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        5);

    qswap(5, QX, 1, QY, 1);

    float128_t* RX = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QX[i].v[0], ==, RX[i].v[0]);
        assert_ullong(QX[i].v[1], ==, RX[i].v[1]);
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qswap_stride(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        9);

    qswap(5, QX, 1, QY, 2);

    float128_t* RX = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QX[i].v[0], ==, RX[i].v[0]);
        assert_ullong(QX[i].v[1], ==, RX[i].v[1]);
    }
    for (uint64_t i = 0; i < 9; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}