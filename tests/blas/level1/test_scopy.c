#include "test.h"

MunitResult test_scopy_all(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);

    scopy(5, (float32_t*)SX, 1, SY, 1, 'n');

    float32_t RY[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){-2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){-4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);

    return MUNIT_OK;
}

MunitResult test_scopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);

    scopy(5, (float32_t*)SX, 2, SY, 1, 'n');

    float32_t RY[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){1.0f}}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);

    return MUNIT_OK;
}

//  N==0 is a no-op: Y is left untouched.
MunitResult test_scopy_zero(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){42.0f}, 1);
    float32_t* SY = svec((float[]){7.0f}, 1);
    scopy(0, SX, 1, SY, 1, 'n');
    float32_t* RY = svec((float[]){7.0f}, 1);
    assert_ulong(SY[0].v, ==, RY[0].v);
    free(SX); free(SY); free(RY);
    return MUNIT_OK;
}
