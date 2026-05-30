#include "test.h"

MunitResult test_saxpy_0(const MunitParameter params[], void *user_data) {
    const uint64_t N = 5;
    const float32_t SA = {*(uint32_t*)&(float){0.0f}};
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);
    float32_t* SY = svec((float[]){10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, 5);

    saxpy(N, SA, SX, 1, SY, 1, 'n');

    float32_t* RY = svec((float[]){10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, 5);

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
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

    saxpy(N, SA, SX, 1, SY, 1, 'n');

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
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

    saxpy(N, SA, SX, 1, SY, 2, 'n');

    for (uint64_t i = 0; i < 9; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
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

    saxpy(N, SA, SX, 1, SY, -1, 'n');

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(SY[i].v, ==, RY[i].v);
    }

    free(SX);
    free(SY);
    free(RY);

    return MUNIT_OK;
}

//  Verifies the rndMode argument actually changes arithmetic. 2^-24 is exactly
//  half a ULP at 1.0, so 1.0 + 2^-24 rounds to 1.0 under round-to-nearest-even
//  (ties to even) but up to the next representable float under round-up. If
//  _set_rounding were a no-op, both modes would give the same answer.
MunitResult test_saxpy_rounding(const MunitParameter params[],
                                void* user_data_or_fixture) {
    const float32_t SA = { SB_REAL32_ONE };
    const float eps = 5.9604644775390625e-8f;  // 2^-24

    //  round toward +inf ('u'): 1.0 + 2^-24 -> next float above 1.0
    float32_t* SXu = svec((float[]){eps}, 1);
    float32_t* SYu = svec((float[]){1.0f}, 1);
    saxpy(1, SA, SXu, 1, SYu, 1, 'u');

    float32_t up = { SB_REAL32_ONE };
    up.v += 1;                              // next representable above 1.0
    assert_ulong(SYu[0].v, ==, up.v);

    //  round to nearest-even ('n'): ties down to 1.0
    float32_t* SXn = svec((float[]){eps}, 1);
    float32_t* SYn = svec((float[]){1.0f}, 1);
    saxpy(1, SA, SXn, 1, SYn, 1, 'n');

    const float32_t one = { SB_REAL32_ONE };
    assert_ulong(SYn[0].v, ==, one.v);

    //  The two modes must disagree, proving the mode took effect.
    assert_ulong(SYu[0].v, !=, SYn[0].v);

    free(SXu); free(SYu); free(SXn); free(SYn);
    return MUNIT_OK;
}
