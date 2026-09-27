#include "test.h"

//  srot: plane rotation with c=0, s=1 maps (x,y) -> (y, -x).
MunitResult test_srot_basic(const MunitParameter params[], void* u) {
    const float32_t c = { 0x00000000 }, s = { 0x3f800000 };  // 0, 1
    float32_t* X = svec((float[]){1.0f, 3.0f}, 2);
    float32_t* Y = svec((float[]){2.0f, 4.0f}, 2);
    srot(2, X, 1, Y, 1, c, s, 'n');
    assert_ulong(X[0].v, ==, 0x40000000u);  // 2
    assert_ulong(X[1].v, ==, 0x40800000u);  // 4
    assert_ulong(Y[0].v, ==, 0xbf800000u);  // -1
    assert_ulong(Y[1].v, ==, 0xc0400000u);  // -3
    free(X); free(Y);
    return MUNIT_OK;
}
//  srotg(a=1, b=0): r=1, z=0, c=1, s=0.
MunitResult test_srotg_basic(const MunitParameter params[], void* u) {
    float32_t a = { 0x3f800000 }, b = { 0x00000000 }, c, s;
    srotg(&a, &b, &c, &s, 'n');
    assert_ulong(a.v, ==, 0x3f800000u);  // r = 1
    assert_ulong(b.v, ==, 0x00000000u);  // z = 0
    assert_ulong(c.v, ==, 0x3f800000u);  // c = 1
    assert_ulong(s.v, ==, 0x00000000u);  // s = 0
    return MUNIT_OK;
}
//  srotg(a=3, b=4): r = 5 (|a|<|b| so roe=b>0 -> positive r).
MunitResult test_srotg_345(const MunitParameter params[], void* u) {
    float32_t a = { 0x40400000 }, b = { 0x40800000 }, c, s;  // 3, 4
    srotg(&a, &b, &c, &s, 'n');
    assert_ulong(a.v, ==, 0x40a00000u);  // r = 5
    return MUNIT_OK;
}
//  srotm flag=0: H = [[1, h12],[h21, 1]] with h21=h12=1. (2,3) -> (5,5).
MunitResult test_srotm_basic(const MunitParameter params[], void* u) {
    float32_t* X = svec((float[]){2.0f}, 1);
    float32_t* Y = svec((float[]){3.0f}, 1);
    float32_t P[5] = {{0x00000000},{0},{0x3f800000},{0x3f800000},{0}};  // flag=0,h21=1,h12=1
    srotm(1, X, 1, Y, 1, P, 'n');
    assert_ulong(X[0].v, ==, 0x40a00000u);  // 5
    assert_ulong(Y[0].v, ==, 0x40a00000u);  // 5
    free(X); free(Y);
    return MUNIT_OK;
}

//  Quad-precision rotation routines. float128_t literals are {lo, hi}; only the
//  hi word carries the exponent/sign for these exact small integers, so tests
//  assert on v[1].
//  qrot c=0, s=1 maps (x,y) -> (y, -x): (2,3) -> (3,-2).
MunitResult test_qrot_basic(const MunitParameter params[], void* u) {
    const float128_t c = {0, 0x0}, s = {0, 0x3fff000000000000};   // 0, 1
    float128_t X = {0, 0x4000000000000000};   // 2
    float128_t Y = {0, 0x4000800000000000};   // 3
    qrot(1, &X, 1, &Y, 1, c, s, 'n');
    assert_ullong(X.v[1], ==, 0x4000800000000000ull);   // 3
    assert_ullong(Y.v[1], ==, 0xc000000000000000ull);   // -2
    return MUNIT_OK;
}
//  qrotg(a=0, b=5): scale=5, roe=b>0 -> r=5 exactly, c=0, s=1.
MunitResult test_qrotg_basic(const MunitParameter params[], void* u) {
    float128_t a = {0, 0x0}, b = {0, 0x4001400000000000}, c, s;   // 0, 5
    qrotg(&a, &b, &c, &s, 'n');
    assert_ullong(a.v[1], ==, 0x4001400000000000ull);   // r = 5
    assert_ullong(c.v[1], ==, 0x0ull);                  // c = 0
    assert_ullong(s.v[1], ==, 0x3fff000000000000ull);   // s = 1
    return MUNIT_OK;
}
//  qrotm flag=0 with h21=h12=1: (2,3) -> (5,5), mirror of test_srotm_basic.
MunitResult test_qrotm_basic(const MunitParameter params[], void* u) {
    float128_t X = {0, 0x4000000000000000};   // 2
    float128_t Y = {0, 0x4000800000000000};   // 3
    float128_t P[5] = {{0,0}, {0,0}, {0, 0x3fff000000000000}, {0, 0x3fff000000000000}, {0,0}};
    qrotm(1, &X, 1, &Y, 1, P, 'n');
    assert_ullong(X.v[1], ==, 0x4001400000000000ull);   // 5
    assert_ullong(Y.v[1], ==, 0x4001400000000000ull);   // 5
    return MUNIT_OK;
}

//  Regression (issue #34): ?rotg must canonicalize its outputs. SoftFloat
//  propagates a NaN input's payload and sign through the arithmetic, so before
//  the fix these came back as e.g. 0x7fef434f (quieted signaling NaN, payload
//  intact) or with the input's sign bit still set. Every NaN the routine
//  writes must now read back as the library's canonical NaN.
//  srotg with a signaling-NaN b: a = 0x2e949262, b = 0x7faf434f.
MunitResult test_srotg_nan_unify(const MunitParameter params[], void* u) {
    float32_t a = { 0x2e949262 }, b = { 0x7faf434f }, c, s;
    srotg(&a, &b, &c, &s, 'n');
    assert_ulong(a.v, ==, (unsigned long)SINGNAN);
    assert_ulong(b.v, ==, (unsigned long)SINGNAN);
    assert_ulong(c.v, ==, (unsigned long)SINGNAN);
    assert_ulong(s.v, ==, (unsigned long)SINGNAN);
    return MUNIT_OK;
}
//  drotg with a signaling NaN: before the fix this came back quieted but with
//  the payload intact (0x7ffc000000000000). (A sign-only NaN like
//  0xfff8000000000000 is not a discriminating input under the 8086
//  specialization this repo builds: f64_abs clears the sign before the NaN
//  reaches the arithmetic, so it canonicalizes even without the unify.)
MunitResult test_drotg_nan_unify(const MunitParameter params[], void* u) {
    float64_t a = { 0x7ff4000000000000 }, b = { 0xbfdc884e0f37fbaa }, c, s;
    drotg(&a, &b, &c, &s, 'n');
    assert_ullong(a.v, ==, (unsigned long long)DOUBNAN);
    assert_ullong(b.v, ==, (unsigned long long)DOUBNAN);
    assert_ullong(c.v, ==, (unsigned long long)DOUBNAN);
    assert_ullong(s.v, ==, (unsigned long long)DOUBNAN);
    return MUNIT_OK;
}
//  hrotg with a signaling NaN (0x7d00), b = 1: quieted to 0x7f00 before the fix.
MunitResult test_hrotg_nan_unify(const MunitParameter params[], void* u) {
    float16_t a = { 0x7d00 }, b = { 0x3c00 }, c, s;
    hrotg(&a, &b, &c, &s, 'n');
    assert_ulong(a.v, ==, (unsigned long)HALFNAN);
    assert_ulong(b.v, ==, (unsigned long)HALFNAN);
    assert_ulong(c.v, ==, (unsigned long)HALFNAN);
    assert_ulong(s.v, ==, (unsigned long)HALFNAN);
    return MUNIT_OK;
}
//  qrotg with a signaling NaN, b = 1: quieted to 0x7fffc000... before the fix.
//  nan_unify_q zeroes the low word as well.
MunitResult test_qrotg_nan_unify(const MunitParameter params[], void* u) {
    float128_t a = {0, 0x7fff400000000000}, b = {0, 0x3fff000000000000}, c, s;
    qrotg(&a, &b, &c, &s, 'n');
    assert_ullong(a.v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(b.v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(c.v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(s.v[1], ==, (unsigned long long)QUADNAN);
    assert_ullong(a.v[0], ==, 0x0ull);
    assert_ullong(b.v[0], ==, 0x0ull);
    assert_ullong(c.v[0], ==, 0x0ull);
    assert_ullong(s.v[0], ==, 0x0ull);
    return MUNIT_OK;
}
