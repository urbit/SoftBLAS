#include "test.h"

//  sdsdot: alpha + sum(x*y) in double precision. 1 + 2*4 + 3*5 = 24.
MunitResult test_sdsdot_basic(const MunitParameter params[], void* u) {
    const float32_t alpha = { SB_REAL32_ONE };
    float32_t* SX = svec((float[]){2.0f, 3.0f}, 2);
    float32_t* SY = svec((float[]){4.0f, 5.0f}, 2);
    float32_t r = sdsdot(2, alpha, SX, 1, SY, 1, 'n');
    assert_ulong(r.v, ==, 0x41c00000u);  // 24
    free(SX); free(SY);
    return MUNIT_OK;
}
//  hsdot: half inputs, single accumulation. 1 + 2*4 + 3*5 = 24.
MunitResult test_hsdot_basic(const MunitParameter params[], void* u) {
    const float16_t alpha = { SB_REAL16_ONE };
    float16_t* X = hvec((uint16_t[]){0x4000, 0x4200}, 2);  // 2, 3
    float16_t* Y = hvec((uint16_t[]){0x4400, 0x4500}, 2);  // 4, 5
    float32_t r = hsdot(2, alpha, X, 1, Y, 1, 'n');
    assert_ulong(r.v, ==, 0x41c00000u);  // 24 (single)
    free(X); free(Y);
    return MUNIT_OK;
}
