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

//  Regression: a NaN or Inf in D1/D2 must NOT hang the rescale loops. The
//  loops test `f_ge(d, gamsq)` = `!(f_lt(d, gamsq))`, which is true for a
//  non-finite d forever (NaN/Inf never compare less-than), so before the fix
//  these spun indefinitely. The loops are now guarded with a finiteness check.
//  These tests passing (completing at all) is the liveness assertion; the
//  value assertions pin the deterministic result.
MunitResult test_srotmg_nan_d1(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x7fc00000 }, d2 = { 0x3f800000 }, x1 = { 0x3f800000 };  // NaN, 1, 1
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5];
    srotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ulong(P[0].v, ==, 0x00000000u);  // flag = 0
    assert_true(nan_test_s(d1));            // NaN propagates to the output
    return MUNIT_OK;
}
MunitResult test_srotmg_inf_d1(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x7f800000 }, d2 = { 0x3f800000 }, x1 = { 0x3f800000 };  // +Inf, 1, 1
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5];
    srotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ulong(d1.v, ==, 0x7f800000u);    // +Inf propagates, no rescale
    return MUNIT_OK;
}
MunitResult test_hrotmg_inf_d1(const MunitParameter params[], void* u) {
    float16_t d1 = { 0x7c00 }, d2 = { 0x3c00 }, x1 = { 0x3c00 };  // +Inf, 1, 1
    const float16_t y1 = { 0x3c00 };
    float16_t P[5];
    hrotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ulong(d1.v, ==, 0x7c00u);        // +Inf propagates, no rescale
    return MUNIT_OK;
}

//  qrotmg(1,1,1,1): quad-precision mirror of test_srotmg_basic.
//  flag = 1; d1,d2 -> 0.5; x1 -> 2. (float128_t is {lo, hi}; hi word checked.)
MunitResult test_qrotmg_basic(const MunitParameter params[], void* u) {
    float128_t d1 = {0, 0x3fff000000000000}, d2 = {0, 0x3fff000000000000},
               x1 = {0, 0x3fff000000000000};
    const float128_t y1 = {0, 0x3fff000000000000};
    float128_t P[5];
    qrotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ullong(d1.v[1],   ==, 0x3ffe000000000000ull);  // 0.5
    assert_ullong(x1.v[1],   ==, 0x4000000000000000ull);  // 2
    assert_ullong(P[0].v[1], ==, 0x3fff000000000000ull);  // flag = 1
    return MUNIT_OK;
}
//  Regression: quad rotmg must not hang on a non-finite D1 (same guard as s/d/h).
MunitResult test_qrotmg_inf_d1(const MunitParameter params[], void* u) {
    float128_t d1 = {0, 0x7fff000000000000}, d2 = {0, 0x3fff000000000000},
               x1 = {0, 0x3fff000000000000};                 // +Inf, 1, 1
    const float128_t y1 = {0, 0x3fff000000000000};
    float128_t P[5];
    qrotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ullong(d1.v[1], ==, 0x7fff000000000000ull);       // +Inf propagates
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
