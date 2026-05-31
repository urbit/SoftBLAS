#include "test.h"

//  C3: N==0 is a no-op across the board. Reductions return canonical 0; the
//  vector ops leave their destination untouched. Single precision is covered in
//  the per-routine files; this sweeps double, half, and quad in one test each.

MunitResult test_dn0(const MunitParameter params[], void* u) {
    const float64_t A = { SB_REAL64_ONE };
    float64_t X[2] = {{ 0x4000000000000000ull }, { 0x4000000000000000ull }};  // 2,2
    float64_t Y[2] = {{ 0x4008000000000000ull }, { 0x4008000000000000ull }};  // 3,3
    //  reductions -> canonical 0
    assert_ullong(dasum(0, X, 1, 'n').v,       ==, (uint64_t)SB_REAL64_ZERO);
    assert_ullong(ddot(0, X, 1, Y, 1, 'n').v,  ==, (uint64_t)SB_REAL64_ZERO);
    assert_ullong(dnrm2(0, X, 1, 'n').v,       ==, (uint64_t)SB_REAL64_ZERO);
    //  vector ops -> destination untouched
    daxpy(0, A, X, 1, Y, 1, 'n');
    dscal(0, A, Y, 1, 'n');
    dcopy(0, X, 1, Y, 1, 'n');
    dswap(0, X, 1, Y, 1, 'n');
    assert_ullong(Y[0].v, ==, 0x4008000000000000ull);
    assert_ullong(Y[1].v, ==, 0x4008000000000000ull);
    assert_ullong(X[0].v, ==, 0x4000000000000000ull);
    return MUNIT_OK;
}

MunitResult test_hn0(const MunitParameter params[], void* u) {
    const float16_t A = { SB_REAL16_ONE };
    float16_t X[2] = {{ 0x4000 }, { 0x4000 }};  // 2,2
    float16_t Y[2] = {{ 0x4200 }, { 0x4200 }};  // 3,3
    assert_ushort(hasum(0, X, 1, 'n').v,      ==, (uint16_t)SB_REAL16_ZERO);
    assert_ushort(hdot(0, X, 1, Y, 1, 'n').v, ==, (uint16_t)SB_REAL16_ZERO);
    assert_ushort(hnrm2(0, X, 1, 'n').v,      ==, (uint16_t)SB_REAL16_ZERO);
    haxpy(0, A, X, 1, Y, 1, 'n');
    hscal(0, A, Y, 1, 'n');
    hcopy(0, X, 1, Y, 1, 'n');
    hswap(0, X, 1, Y, 1, 'n');
    assert_ushort(Y[0].v, ==, 0x4200);
    assert_ushort(X[0].v, ==, 0x4000);
    return MUNIT_OK;
}

MunitResult test_qn0(const MunitParameter params[], void* u) {
    const float128_t A = {{ 0x0ull, 0x3fff000000000000ull }};  // 1.0
    float128_t X[2] = {{{ 0x0ull, 0x4000000000000000ull }},    // 2,2
                       {{ 0x0ull, 0x4000000000000000ull }}};
    float128_t Y[2] = {{{ 0x0ull, 0x4000800000000000ull }},    // 3,3
                       {{ 0x0ull, 0x4000800000000000ull }}};
    float128_t r;
    r = qasum(0, X, 1, 'n');       assert_ullong(r.v[1], ==, 0x0ull); assert_ullong(r.v[0], ==, 0x0ull);
    r = qdot(0, X, 1, Y, 1, 'n');  assert_ullong(r.v[1], ==, 0x0ull); assert_ullong(r.v[0], ==, 0x0ull);
    r = qnrm2(0, X, 1, 'n');       assert_ullong(r.v[1], ==, 0x0ull); assert_ullong(r.v[0], ==, 0x0ull);
    qaxpy(0, A, X, 1, Y, 1, 'n');
    qscal(0, A, Y, 1, 'n');
    qcopy(0, X, 1, Y, 1, 'n');
    qswap(0, X, 1, Y, 1, 'n');
    assert_ullong(Y[0].v[1], ==, 0x4000800000000000ull);
    assert_ullong(X[0].v[1], ==, 0x4000000000000000ull);
    return MUNIT_OK;
}

//  B5: a stride large enough that (N-1)*incX overflows uint64_t must be
//  rejected up front (return canonical 0) rather than walking off the buffer.
//  The single-element buffer would be a wild OOB read without the guard; this
//  test passing under ASan confirms nothing is dereferenced.
MunitResult test_nrm2_huge_stride(const MunitParameter params[], void* u) {
    float32_t X[1] = {{ 0x3f800000u }};   // 1 element only
    const uint64_t huge = 0x1000000000ull;  // 2^36; (huge-1)*huge ~ 2^72 wraps
    assert_ulong(snrm2(huge, X, huge, 'n').v,       ==, (uint32_t)SB_REAL32_ZERO);
    assert_ulong(snrm2_B(huge, X, huge, 'n').v,     ==, (uint32_t)SB_REAL32_ZERO);
    assert_ulong(sasum(huge, X, huge, 'n').v,       ==, (uint32_t)SB_REAL32_ZERO);
    return MUNIT_OK;
}
