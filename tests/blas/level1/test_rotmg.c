#include "test.h"

//  srotmg(d1=d2=x1=y1=1): the no-rescale flag=1 branch.
//  u = 1 + h11*h22 = 2; d1,d2 -> 0.5; x1 -> y1*u = 2; h11 = h22 = 1.
MunitResult test_srotmg_basic(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x3f800000 }, d2 = { 0x3f800000 }, x1 = { 0x3f800000 };
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5];
    srotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ulong(d1.v,   ==, 0x3f000000u);  // 0.5
    assert_ulong(d2.v,   ==, 0x3f000000u);  // 0.5
    assert_ulong(x1.v,   ==, 0x40000000u);  // 2
    assert_ulong(P[0].v, ==, 0x3f800000u);  // flag = 1
    assert_ulong(P[1].v, ==, 0x3f800000u);  // h11 = 1
    assert_ulong(P[4].v, ==, 0x3f800000u);  // h22 = 1
    return MUNIT_OK;
}
//  srotmg with d2*y1 == 0 returns the "no rotation" flag (-2) immediately.
MunitResult test_srotmg_flag_neg2(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x3f800000 }, d2 = { 0x00000000 }, x1 = { 0x3f800000 };
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5] = {{0},{0},{0},{0},{0}};
    srotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ulong(P[0].v, ==, 0xc0000000u);  // flag = -2
    return MUNIT_OK;
}
//  drotmg(1,1,1,1): same as srotmg but double precision.
MunitResult test_drotmg_basic(const MunitParameter params[], void* u) {
    float64_t d1 = { 0x3ff0000000000000 }, d2 = { 0x3ff0000000000000 }, x1 = { 0x3ff0000000000000 };
    const float64_t y1 = { 0x3ff0000000000000 };
    float64_t P[5];
    drotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ullong(d1.v,   ==, 0x3fe0000000000000ull);  // 0.5
    assert_ullong(x1.v,   ==, 0x4000000000000000ull);  // 2
    assert_ullong(P[0].v, ==, 0x3ff0000000000000ull);  // flag = 1
    return MUNIT_OK;
}

//  Singular case (q2 < 0): zeroes the result and sets flag = -1. This is the
//  path that previously used `goto store`; it must still produce all zeros.
MunitResult test_srotmg_singular(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x3f800000 }, d2 = { 0xbf800000 }, x1 = { 0x3f800000 };  // 1, -1, 1
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5] = {{0},{0},{0},{0},{0}};
    srotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ulong(P[0].v, ==, 0xbf800000u);  // flag = -1
    assert_ulong(d1.v, ==, 0x00000000u);
    assert_ulong(d2.v, ==, 0x00000000u);
    assert_ulong(x1.v, ==, 0x00000000u);
    for (int i = 1; i <= 4; i++) assert_ulong(P[i].v, ==, 0x00000000u);
    return MUNIT_OK;
}
