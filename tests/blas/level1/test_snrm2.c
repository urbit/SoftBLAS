#include "test.h"

MunitResult test_snrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){0.0f, 0.0f}, 2);

    float32_t S = (float32_t) snrm2(2, (float32_t*)SX, 1);
    float32_t R = (float32_t){0.0f};

    assert_ulong(S.v, ==, R.v);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_snrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, -2.0f, 3.0f, -4.0f, 5.0f}, 5);

    float32_t S = (float32_t) snrm2(5, (float32_t*)SX, 1);
    float32_t R = {*(uint32_t*)&(float){7.416198487f}};  // sqrt(55.0)

    assert_ulong(S.v, ==, R.v);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_snrm2_stride(const MunitParameter params[],
                              void *user_data) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f}, 9);

    float32_t S = (float32_t) snrm2(5, (float32_t*)SX, 2);
    float32_t R = {*(uint32_t*)&(float){2.236067977f}};  // sqrt(5.0)

    assert_ulong(S.v, ==, R.v);

    free(SX);

    return MUNIT_OK;
}
