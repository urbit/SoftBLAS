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
