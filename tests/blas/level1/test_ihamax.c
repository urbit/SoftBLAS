#include "test.h"

MunitResult test_ihamax_0(const MunitParameter params[],
                          void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x0 , 0x0 , 0x0 }, 3);

    uint64_t I = ihamax(3, HX, 1);
    uint64_t R = 0;

    assert_int(I, ==, R);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_ihamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);

    uint64_t I = ihamax(5, HX, 1);
    uint64_t R = 4;

    assert_int(I, ==, R);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_ihamax_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);

    uint64_t I = ihamax(2, HX, 3);
    // note that ihamax counts against the stride index, not the overall index
    uint64_t R = 1;

    assert_int(I, ==, R);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_ihamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4200 , 0x4500 , 0x4400 , 0x4000 }, 5);

    uint64_t I = ihamax(5, HX, 1);
    uint64_t R = 2;

    assert_int(I, ==, R);
    
    free(HX);

    return MUNIT_OK;
}
