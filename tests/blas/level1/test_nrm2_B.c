#include "test.h"

//  Anderson safe-scaling nrm2: correct in the medium band, and -- the whole
//  point -- overflow-safe for large magnitudes where naive sqrt(sum x^2) would
//  overflow.
MunitResult test_snrm2_B_med(const MunitParameter params[], void* u) {
    float32_t* X = svec((float[]){3.0f, 4.0f}, 2);
    assert_ulong(snrm2_B(2, X, 1, 'n').v, ==, 0x40a00000u);   // 5
    free(X); return MUNIT_OK;
}
MunitResult test_snrm2_B_big(const MunitParameter params[], void* u) {
    float32_t X[1] = {{ 0x7b800000 }};                        // 2^120 (2^240 overflows)
    assert_ulong(snrm2_B(1, X, 1, 'n').v, ==, 0x7b800000u);   // 2^120
    return MUNIT_OK;
}
MunitResult test_snrm2_B_small(const MunitParameter params[], void* u) {
    float32_t X[1] = {{ 0x0d800000 }};                        // 2^-100 (2^-200 underflows)
    assert_ulong(snrm2_B(1, X, 1, 'n').v, ==, 0x0d800000u);   // 2^-100
    return MUNIT_OK;
}
MunitResult test_dnrm2_B_med(const MunitParameter params[], void* u) {
    float64_t* X = dvec((double[]){3.0, 4.0}, 2);
    assert_ullong(dnrm2_B(2, X, 1, 'n').v, ==, 0x4014000000000000ull);  // 5
    free(X); return MUNIT_OK;
}
MunitResult test_dnrm2_B_big(const MunitParameter params[], void* u) {
    float64_t X[1] = {{ 0x7e70000000000000ull }};             // 2^1000 (2^2000 overflows)
    assert_ullong(dnrm2_B(1, X, 1, 'n').v, ==, 0x7e70000000000000ull);
    return MUNIT_OK;
}
MunitResult test_hnrm2_B_med(const MunitParameter params[], void* u) {
    float16_t X[2] = {{ 0x4200 }, { 0x4400 }};                // 3, 4 -> 5
    assert_ushort(hnrm2_B(2, X, 1, 'n').v, ==, 0x4500);       // 5
    return MUNIT_OK;
}
MunitResult test_hnrm2_B_big(const MunitParameter params[], void* u) {
    float16_t X[1] = {{ 0x6400 }};                            // 1024 (1024^2 overflows float16)
    assert_ushort(hnrm2_B(1, X, 1, 'n').v, ==, 0x6400);       // 1024
    return MUNIT_OK;
}
MunitResult test_qnrm2_B_med(const MunitParameter params[], void* u) {
    float128_t* X = qvec((float128_pair_t[]){
        {.hi=0x4000800000000000,.lo=0},{.hi=0x4001000000000000,.lo=0}}, 2);  // 3,4
    float128_t r = qnrm2_B(2, X, 1, 'n');
    assert_ullong(r.v[1], ==, 0x4001400000000000ull);  // 5
    assert_ullong(r.v[0], ==, 0ull);
    free(X); return MUNIT_OK;
}
