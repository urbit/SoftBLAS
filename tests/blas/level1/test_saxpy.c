#include "test.h"

//  Test SAXPY, A*X+Y.
//
//  void saxpy(uint64_t N, float32_t SA, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
//
MunitResult test_saxpy_0(const MunitParameter params[],
                         void *user_data) {
    const float32_t SA = { 0.0 };
    union sing SX[] = {
        {.s = {*(float*)&(float){1.0f}}},
        {.s = {*(float*)&(float){-2.0f}}},
        {.s = {*(float*)&(float){3.0f}}},
        {.s = {*(float*)&(float){-4.0f}}},
        {.s = {*(float*)&(float){5.0f}}}
    };
    union sing SY[] = {
        {.s = {*(float*)&(float){10.0f}}},
        {.s = {*(float*)&(float){20.0f}}},
        {.s = {*(float*)&(float){30.0f}}},
        {.s = {*(float*)&(float){40.0f}}},
        {.s = {*(float*)&(float){50.0f}}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    union sing RY[] = {
        {.s = {*(float*)&(float){10.0f}}},
        {.s = {*(float*)&(float){20.0f}}},
        {.s = {*(float*)&(float){30.0f}}},
        {.s = {*(float*)&(float){40.0f}}},
        {.s = {*(float*)&(float){50.0f}}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].u, ==, RY[i].u);
    }

    return MUNIT_OK;
}

MunitResult test_saxpy_sum(const MunitParameter params[],
                           void *user_data) {
    const float32_t SA = { 2.0 };
    union sing SX[] = {
        {.s = {*(float*)&(float){1.0f}}},
        {.s = {*(float*)&(float){-2.0f}}},
        {.s = {*(float*)&(float){3.0f}}},
        {.s = {*(float*)&(float){-4.0f}}},
        {.s = {*(float*)&(float){5.0f}}}
    };
    union sing SY[] = {
        {.s = {*(float*)&(float){10.0f}}},
        {.s = {*(float*)&(float){20.0f}}},
        {.s = {*(float*)&(float){30.0f}}},
        {.s = {*(float*)&(float){40.0f}}},
        {.s = {*(float*)&(float){50.0f}}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    union sing RY[] = {
        {.s = {*(float*)&(float){12.0f}}},
        {.s = {*(float*)&(float){16.0f}}},
        {.s = {*(float*)&(float){24.0f}}},
        {.s = {*(float*)&(float){32.0f}}},
        {.s = {*(float*)&(float){60.0f}}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].u, ==, RY[i].u);
    }

    return MUNIT_OK;
}
