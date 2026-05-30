#include "test.h"

MunitResult test_sscal_0(const MunitParameter params[], void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){0.0f}};
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);

    sscal(5, SA, SX, 1, 'n');

    float32_t* RX = svec((float[]){0.0f, -0.0f, 0.0f, -0.0f, 0.0f}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SX[i].v, ==, RX[i].v);
    }

    free(SX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_sscal_12345(const MunitParameter params[],
                           void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){10.5f}};
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);

    sscal(5, SA, SX, 1, 'n');

    float32_t* RX = svec((float[]){10.5f, -21.0f, 31.5f, -42.0f, 52.5f}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SX[i].v, ==, RX[i].v);
    }

    free(SX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_sscal_stride(const MunitParameter params[],
                              void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){2.5f}};
    float32_t* SX = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);

    sscal(5, SA, SX, 2, 'n');

    float32_t* RX = svec((float[]){2.5f, 2.0f, 2.5f, 2.0f, 2.5f, 2.0f, 2.5f, 2.0f, 2.5f}, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ulong(SX[i].v, ==, RX[i].v);
    }

    free(SX);
    free(RX);

    return MUNIT_OK;
}

//  N==0 is a no-op: X is left untouched.
MunitResult test_sscal_zero(const MunitParameter params[],
                            void* user_data_or_fixture) {
    const float32_t SA = { SB_REAL32_ZERO };
    float32_t* SX = svec((float[]){42.0f}, 1);
    sscal(0, SA, SX, 1, 'n');
    float32_t* RX = svec((float[]){42.0f}, 1);
    assert_ulong(SX[0].v, ==, RX[0].v);
    free(SX); free(RX);
    return MUNIT_OK;
}
