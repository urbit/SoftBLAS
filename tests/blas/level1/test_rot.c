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
