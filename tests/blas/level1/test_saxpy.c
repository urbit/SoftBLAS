#include "test.h"

MunitResult test_saxpy_0(const MunitParameter params[], void *user_data) {
    const uint64_t N = 5;
    const float32_t SA = {*(uint32_t*)&(float){0.0f}};
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, 5);

    saxpy(N, SA, SX, 1, SY, 1);

    float32_t* RY = svec((float[]){10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, 5);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_saxpy_sum(const MunitParameter params[],
                           void *user_data) {
    const uint64_t N = 5;
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, N);
    float32_t* SY = svec((float[]){10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, N);
    float32_t* RY = svec((float[]){12.0f, 16.0f, 36.0f, 32.0f, 60.0f}, N);

    saxpy(N, SA, SX, 1, SY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_saxpy_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, N);
    float32_t* SY = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);
    float32_t* RY = svec((float[]){3.0f, 2.0f, 5.0f, 2.0f, 7.0f, 2.0f, 9.0f, 2.0f, 11.0f}, 9);

    saxpy(N, SA, SX, 1, SY, 2);

    for (uint64_t i = 0; i < 9; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_saxpy_neg_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, N);
    float32_t* SY = svec((float[]){5.0f, 4.0f, 3.0f, 2.0f, 1.0f}, N);
    float32_t* RY = svec((float[]){15.0f, 12.0f, 9.0f, 6.0f, 3.0f}, N);

    saxpy(N, SA, SX, 1, SY, -1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}
