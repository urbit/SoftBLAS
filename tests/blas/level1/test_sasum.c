#include "test.h"

MunitResult test_sasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    const uint64_t N = 2;
    float32_t* SX = svec((float[]){0.0f, 0.0f}, N);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1);
    float32_t R = (float32_t){0.0f};

    assert_int(S.v, ==, R.v);
    
    return MUNIT_OK;
}

MunitResult test_sasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    const uint64_t N = 5;
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, N);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1);
    float32_t R = {*(uint32_t*)&(float){15.0f}};
    
    assert_int(S.v, ==, R.v);
    
    return MUNIT_OK;
}

// TODO test stride
