#include "test.h"

MunitResult test_sdot_0(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);

    float32_t SZ = sdot(5, (float32_t*)SX, 1, SY, 1, 'n');
    float32_t RZ = {*(uint32_t*)&(float){0.0f}};

    assert_ulong(SZ.v, ==, RZ.v);

    free(SX);
    free(SY);

    return MUNIT_OK;
}

MunitResult test_sdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){5.0f, 4.0f, 3.0f, 2.0f, 1.0f}, 5);

    float32_t SZ = sdot(5, (float32_t*)SX, 1, SY, 1, 'n');
    float32_t RZ = {*(uint32_t*)&(float){35.0f}};

    assert_ulong(SZ.v, ==, RZ.v);

    free(SX);
    free(SY);

    return MUNIT_OK;
}

MunitResult test_sdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);

    float32_t SZ = sdot(5, (float32_t*)SX, 1, SY, 2, 'n');
    float32_t RZ = {*(uint32_t*)&(float){15.0f}};

    assert_ulong(SZ.v, ==, RZ.v);

    free(SX);
    free(SY);

    return MUNIT_OK;
}

MunitResult test_sdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){5.0f, 4.0f, 3.0f, 2.0f, 1.0f}, 5);

    float32_t SZ = sdot(5, (float32_t*)SX, 1, SY, -1, 'n');
    float32_t RZ = {*(uint32_t*)&(float){55.0f}};

    assert_ulong(SZ.v, ==, RZ.v);

    free(SX);
    free(SY);

    return MUNIT_OK;
}

//  N==0 is a no-op: returns 0.
MunitResult test_sdot_zero(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){42.0f}, 1);
    float32_t* SY = svec((float[]){7.0f}, 1);
    float32_t r = sdot(0, SX, 1, SY, 1, 'n');
    assert_ulong(r.v, ==, (uint32_t)SB_REAL32_ZERO);
    free(SX); free(SY);
    return MUNIT_OK;
}

//  C6: all five rounding modes through the dot accumulation. dot([1, 2^-24],
//  [1, 1]) = 1 + 2^-24, an exact tie; the negative case distinguishes z/d/u/a.
MunitResult test_sdot_rounding_modes(const MunitParameter params[], void* u) {
    const uint32_t ONE = 0x3f800000u, UP = 0x3f800001u,
                   NEG_ONE = 0xbf800000u, NEG_UP = 0xbf800001u;
    struct { char m; uint32_t want; } pos[] = {
        {'n', ONE}, {'z', ONE}, {'d', ONE}, {'u', UP}, {'a', UP} };
    for (uint64_t k = 0; k < 5; k++) {
        float32_t* X = svec((float[]){ 1.0f, 0x1p-24f }, 2);
        float32_t* Y = svec((float[]){ 1.0f, 1.0f }, 2);
        assert_ulong(sdot(2, X, 1, Y, 1, pos[k].m).v, ==, pos[k].want);
        free(X); free(Y);
    }
    struct { char m; uint32_t want; } neg[] = {
        {'n', NEG_ONE}, {'z', NEG_ONE}, {'u', NEG_ONE}, {'d', NEG_UP}, {'a', NEG_UP} };
    for (uint64_t k = 0; k < 5; k++) {
        float32_t* X = svec((float[]){ -1.0f, -0x1p-24f }, 2);
        float32_t* Y = svec((float[]){ 1.0f, 1.0f }, 2);
        assert_ulong(sdot(2, X, 1, Y, 1, neg[k].m).v, ==, neg[k].want);
        free(X); free(Y);
    }
    return MUNIT_OK;
}
