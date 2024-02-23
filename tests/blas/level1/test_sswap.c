#include "test.h"

MunitResult test_sswap_two(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);

    sswap(5, SX, 1, SY, 1);

    float32_t* RX = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);
    float32_t* RY = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SX[i].v, ==, RX[i].v);
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_sswap_stride(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 9);

    sswap(5, SX, 1, SY, 2);

    float32_t* RX = svec((float[]){0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, 5);
    float32_t* RY = svec((float[]){1.0f, 0.0f, -2.0f, 0.0f, 3.0f, 0.0f, -4.0f, 0.0f, 5.0f}, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SX[i].v, ==, RX[i].v);
    }
    for (uint64_t i = 0; i < 9; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}
