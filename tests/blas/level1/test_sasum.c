#include "test.h"
#include <stdio.h>

//  Test SASUM, sum of absolute values.
//
//  float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX);
//
MunitResult test_sasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float32_t SX[] = {
        {*(uint32_t*)&(float){0.0f}},
        {*(uint32_t*)&(float){0.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1);
    float32_t R = (float32_t){0.0f};

    assert_int(S.v, ==, R.v);
    
    return MUNIT_OK;
}

MunitResult test_sasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float32_t SX[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){-2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){-4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);

    float32_t S = (float32_t) sasum(N, (float32_t*)SX, 1);
    float32_t R = {*(uint32_t*)&(float){15.0f}};
    
    assert_int(S.v, ==, R.v);
    
    return MUNIT_OK;
}

// TODO test stride
