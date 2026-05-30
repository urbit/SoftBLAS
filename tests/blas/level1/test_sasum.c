#include "test.h"

MunitResult test_sasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    const uint64_t N = 2;
    float32_t* SX = svec((float[]){0.0f, 0.0f}, N);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1, 'n');
    float32_t R = (float32_t){0.0f};

    assert_ulong(S.v, ==, R.v);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_sasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    const uint64_t N = 5;
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, N);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1, 'n');
    float32_t R = {*(uint32_t*)&(float){15.0f}};
    
    assert_ulong(S.v, ==, R.v);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_sasum_stride(const MunitParameter params[],
                              void *user_data) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);

    float32_t S = (float32_t) sasum(5, (float32_t*)SX, 2, 'n');
    float32_t R = {*(uint32_t*)&(float){5.0f}};

    assert_ulong(S.v, ==, R.v);

    free(SX);

    return MUNIT_OK;
}

//  N==0 is a no-op: returns 0, reads nothing.
MunitResult test_sasum_zero(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){42.0f}, 1);
    float32_t r = sasum(0, SX, 1, 'n');
    assert_ulong(r.v, ==, (uint32_t)SB_REAL32_ZERO);
    free(SX);
    return MUNIT_OK;
}
