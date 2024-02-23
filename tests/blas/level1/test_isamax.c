#include "test.h"

MunitResult test_isamax_0(const MunitParameter params[],
                          void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){0.0f, 0.0f, 0.0f}, 3);

    uint64_t I = isamax(3, SX, 1);
    uint64_t R = 0;

    assert_int(I, ==, R);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_isamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);

    uint64_t I = isamax(5, SX, 1);
    uint64_t R = 4;

    assert_int(I, ==, R);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_isamax_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, 5);

    uint64_t I = isamax(2, SX, 3);
    // note that isamax counts against the stride index, not the overall index
    uint64_t R = 1;

    assert_int(I, ==, R);
    
    free(SX);

    return MUNIT_OK;
}

MunitResult test_isamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float32_t* SX = svec((float[]){1.0f, 3.0f, 5.0f, 4.0f, 2.0f}, 5);

    uint64_t I = isamax(5, SX, 1);
    uint64_t R = 2;

    assert_int(I, ==, R);
    
    free(SX);

    return MUNIT_OK;
}
