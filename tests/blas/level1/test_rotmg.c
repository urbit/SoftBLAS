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

//  Regression (issue #33): a negative D1 must not hang the first rescale loop.
//  d2 < 0 with a y1 small enough that q2 = (d2*y1)*y1 underflows to -0 slips
//  past the `q2 < 0` singularity test, and the flag = 1 branch then swaps that
//  negative d2 into d1. The old `d1 != 0 && d1 <= rgamsq` guard never went
//  false for a negative d1 -- multiplying by gamsq keeps it negative all the
//  way to -inf, which is still <= rgamsq -- so the loop spun forever. The
//  guard is now `0 < d1`. Completing at all is the liveness assertion; the
//  value assertions pin the (unrescaled) flag = 1 result.
MunitResult test_srotmg_neg_d1_underflow(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x3f800000 }, d2 = { 0xbf800000 }, x1 = { 0x00000000 };  // 1, -1, 0
    const float32_t y1 = { 0x0d800000 };                                      // 2^-100
    float32_t P[5] = {{0},{0},{0},{0},{0}};
    srotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ulong(P[0].v, ==, 0x3f800000u);  // flag = 1
    assert_ulong(d1.v,   ==, 0xbf800000u);  // d1 = d2/u = -1 (swapped in)
    assert_ulong(d2.v,   ==, 0x3f800000u);  // d2 = d1/u = 1
    assert_ulong(x1.v,   ==, 0x0d800000u);  // x1 = y1*u = 2^-100
    assert_ulong(P[1].v, ==, 0x80000000u);  // h11 = p1/p2 = 0/-2^-100 = -0
    assert_ulong(P[4].v, ==, 0x00000000u);  // h22 = x1/y1 = 0
    return MUNIT_OK;
}
MunitResult test_drotmg_neg_d1_underflow(const MunitParameter params[], void* u) {
    float64_t d1 = { 0x3ff0000000000000 }, d2 = { 0xbff0000000000000 },
              x1 = { 0x0000000000000000 };                      // 1, -1, 0
    const float64_t y1 = { 0x1a70000000000000 };                // 2^-600
    float64_t P[5] = {{0},{0},{0},{0},{0}};
    drotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ullong(P[0].v, ==, 0x3ff0000000000000ull);  // flag = 1
    assert_ullong(d1.v,   ==, 0xbff0000000000000ull);  // -1
    assert_ullong(d2.v,   ==, 0x3ff0000000000000ull);  //  1
    assert_ullong(x1.v,   ==, 0x1a70000000000000ull);  // 2^-600
    assert_ullong(P[1].v, ==, 0x8000000000000000ull);  // -0
    assert_ullong(P[4].v, ==, 0x0000000000000000ull);
    return MUNIT_OK;
}
MunitResult test_hrotmg_neg_d1_underflow(const MunitParameter params[], void* u) {
    float16_t d1 = { 0x3c00 }, d2 = { 0xbc00 }, x1 = { 0x0000 };  // 1, -1, 0
    const float16_t y1 = { 0x0800 };                              // 2^-13
    float16_t P[5] = {{0},{0},{0},{0},{0}};
    hrotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ulong(P[0].v, ==, 0x3c00u);   // flag = 1
    assert_ulong(d1.v,   ==, 0xbc00u);   // -1
    assert_ulong(d2.v,   ==, 0x3c00u);   //  1
    assert_ulong(x1.v,   ==, 0x0800u);   // 2^-13
    assert_ulong(P[1].v, ==, 0x8000u);   // -0
    assert_ulong(P[4].v, ==, 0x0000u);
    return MUNIT_OK;
}
MunitResult test_qrotmg_neg_d1_underflow(const MunitParameter params[], void* u) {
    float128_t d1 = {0, 0x3fff000000000000}, d2 = {0, 0xbfff000000000000},
               x1 = {0, 0x0000000000000000};                     // 1, -1, 0
    const float128_t y1 = {0, 0x1cd7000000000000};               // 2^-9000
    float128_t P[5] = {{0,0},{0,0},{0,0},{0,0},{0,0}};
    qrotmg(&d1, &d2, &x1, y1, P, 'n');   // must return, not hang
    assert_ullong(P[0].v[1], ==, 0x3fff000000000000ull);  // flag = 1
    assert_ullong(d1.v[1],   ==, 0xbfff000000000000ull);  // -1
    assert_ullong(d2.v[1],   ==, 0x3fff000000000000ull);  //  1
    assert_ullong(x1.v[1],   ==, 0x1cd7000000000000ull);  // 2^-9000
    assert_ullong(P[1].v[1], ==, 0x8000000000000000ull);  // -0
    assert_ullong(P[4].v[1], ==, 0x0000000000000000ull);
    return MUNIT_OK;
}

//  Regression (issue #34): ?rotmg must canonicalize the values it writes. A
//  NaN x1 drives NaNs into D1/D2/X1 and into the h21/h12 it stores for
//  flag = 0; each must read back as the canonical NaN rather than carrying
//  SoftFloat's propagated payload and sign. (P[1] and P[4] are not written on
//  the flag = 0 path, so they keep the caller's incoming values.)
MunitResult test_srotmg_nan_unify(const MunitParameter params[], void* u) {
    float32_t d1 = { 0x3f800000 }, d2 = { 0x3f800000 }, x1 = { 0xffc00001 };  // 1, 1, -NaN
    const float32_t y1 = { 0x3f800000 };
    float32_t P[5] = {{0},{0},{0},{0},{0}};
    srotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ulong(P[0].v, ==, 0x00000000u);            // flag = 0
    assert_ulong(d1.v,   ==, (unsigned long)SINGNAN);
    assert_ulong(d2.v,   ==, (unsigned long)SINGNAN);
    assert_ulong(x1.v,   ==, (unsigned long)SINGNAN);
    assert_ulong(P[2].v, ==, (unsigned long)SINGNAN);  // h21
    assert_ulong(P[3].v, ==, (unsigned long)SINGNAN);  // h12
    return MUNIT_OK;
}
MunitResult test_drotmg_nan_unify(const MunitParameter params[], void* u) {
    float64_t d1 = { 0x3ff0000000000000 }, d2 = { 0x3ff0000000000000 },
              x1 = { 0xfff8000000000001 };                        // 1, 1, -NaN
    const float64_t y1 = { 0x3ff0000000000000 };
    float64_t P[5] = {{0},{0},{0},{0},{0}};
    drotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ullong(P[0].v, ==, 0x0ull);                      // flag = 0
    assert_ullong(d1.v,   ==, (unsigned long long)DOUBNAN);
    assert_ullong(d2.v,   ==, (unsigned long long)DOUBNAN);
    assert_ullong(x1.v,   ==, (unsigned long long)DOUBNAN);
    assert_ullong(P[2].v, ==, (unsigned long long)DOUBNAN);
    assert_ullong(P[3].v, ==, (unsigned long long)DOUBNAN);
    return MUNIT_OK;
}
MunitResult test_hrotmg_nan_unify(const MunitParameter params[], void* u) {
    float16_t d1 = { 0x3c00 }, d2 = { 0x3c00 }, x1 = { 0xfe01 };  // 1, 1, -NaN
    const float16_t y1 = { 0x3c00 };
    float16_t P[5] = {{0},{0},{0},{0},{0}};
    hrotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ulong(P[0].v, ==, 0x0000u);                // flag = 0
    assert_ulong(d1.v,   ==, (unsigned long)HALFNAN);
    assert_ulong(d2.v,   ==, (unsigned long)HALFNAN);
    assert_ulong(x1.v,   ==, (unsigned long)HALFNAN);
    assert_ulong(P[2].v, ==, (unsigned long)HALFNAN);
    assert_ulong(P[3].v, ==, (unsigned long)HALFNAN);
    return MUNIT_OK;
}
MunitResult test_qrotmg_nan_unify(const MunitParameter params[], void* u) {
    float128_t d1 = {0, 0x3fff000000000000}, d2 = {0, 0x3fff000000000000},
               x1 = {1, 0xffff800000000000};                      // 1, 1, -NaN
    const float128_t y1 = {0, 0x3fff000000000000};
    float128_t P[5] = {{0,0},{0,0},{0,0},{0,0},{0,0}};
    qrotmg(&d1, &d2, &x1, y1, P, 'n');
    assert_ullong(P[0].v[1], ==, 0x0ull);                      // flag = 0
    assert_ullong(d1.v[1],   ==, (unsigned long long)QUADNAN);
    assert_ullong(d2.v[1],   ==, (unsigned long long)QUADNAN);
    assert_ullong(x1.v[1],   ==, (unsigned long long)QUADNAN);
    assert_ullong(P[2].v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(P[3].v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(d1.v[0],   ==, 0x0ull);   // nan_unify_q clears the low word
    assert_ullong(x1.v[0],   ==, 0x0ull);
    return MUNIT_OK;
}
