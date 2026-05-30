#include "test.h"

MunitResult test_qaxpy_0(const MunitParameter params[], void *user_data) {
    // XXX note that 0,1 is lo,hi
    const float128_t QA = {0x0000000000000000, 0x0000000000000000};
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4002400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003e00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004900000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4002400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003e00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004900000000000, .lo = 0x0000000000000000}},
        5);

    qaxpy(5, QA, QX, 1, QY, 1, 'n');

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qaxpy_sum(const MunitParameter params[], void *user_data) {
    // XXX note that 0,1 is lo,hi
    const float128_t QA = {0x0000000000000000, 0x4000000000000000};
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},
            {.hi = 0xc000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0xc001000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x4002400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003e00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004900000000000, .lo = 0x0000000000000000}},
        5);
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4002800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4003000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004200000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4004e00000000000, .lo = 0x0000000000000000}},
        5);

    qaxpy(5, QA, QX, 1, QY, 1, 'n');

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qaxpy_stride(const MunitParameter params[], void *user_data) {
    // XXX note that 0,1 is lo,hi
    const float128_t QA = {0x0000000000000000, 0x4000000000000000};
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
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001400000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001c00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002200000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002600000000000, .lo = 0x0000000000000000}},
        9);

    qaxpy(5, QA, QX, 1, QY, 2, 'n');

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_qaxpy_neg_stride(const MunitParameter params[], void *user_data) {
    // XXX note that 0,1 is lo,hi
    const float128_t QA = {0x0000000000000000, 0x4000000000000000};
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
    float128_t* RY = qvec((float128_pair_t[]){
            {.hi = 0x4002e00000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4002200000000000, .lo = 0x0000000000000000},
            {.hi = 0x4001800000000000, .lo = 0x0000000000000000},
            {.hi = 0x4000800000000000, .lo = 0x0000000000000000}},
        5);

    qaxpy(5, QA, QX, 1, QY, -1, 'n');

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(QY[i].v[0], ==, RY[i].v[0]);
        assert_ullong(QY[i].v[1], ==, RY[i].v[1]);
    }

    free(QX);
    free(QY);
    free(RY);

    return MUNIT_OK;
}



//  qaxpy must canonicalize NaNs in EVERY written element, not just QY[0].
//  QX[1] is a non-canonical NaN whose (quieted) payload SoftFloat propagates;
//  without a per-element nan_unify_q, QY[1] escapes non-canonical and breaks
//  the determinism guarantee.
MunitResult test_qaxpy_nan_unify(const MunitParameter params[],
                                 void* user_data_or_fixture) {
    float128_t QA = {{ 0x0000000000000000, 0x3fff000000000000 }};  // 1.0
    float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3fff000000000000, .lo = 0x0000000000000000},   // 1.0
            {.hi = 0x7fff000000000001, .lo = 0x00000000deadbeef}},  // non-canonical NaN
        2);
    float128_t* QY = qvec((float128_pair_t[]){
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000},
            {.hi = 0x0000000000000000, .lo = 0x0000000000000000}},
        2);

    qaxpy(2, QA, QX, 1, QY, 1, 'n');

    //  QY[1] = 1*NaN + 0 -> canonical QUADNAN with zeroed low word.
    assert_ullong(QY[1].v[1], ==, QUADNAN);
    assert_ullong(QY[1].v[0], ==, 0x0000000000000000);

    free(QX); free(QY);
    return MUNIT_OK;
}
