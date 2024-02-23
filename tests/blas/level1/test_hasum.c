#include "test.h"

MunitResult test_hasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    const uint64_t N = 2;
    float16_t* HX = hvec((uint16_t[]){0x0 , 0x0 }, N);

    float16_t H = (float16_t) hasum(N, (float16_t*)HX, 1);
    float16_t R = (float16_t){0x0 };

    assert_ushort(H.v, ==, R.v);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_hasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    const uint64_t N = 5;
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, N);

    float16_t H = (float16_t) hasum(N, (float16_t*)HX, 1);
    float16_t R = {*(uint16_t*)&(uint16_t){0x4b80 }};
    
    assert_ushort(H.v, ==, R.v);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_hasum_stride(const MunitParameter params[],
                              void *user_data) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 }, 9);

    float16_t H = (float16_t) hasum(5, (float16_t*)HX, 2);
    float16_t R = {*(uint16_t*)&(uint16_t){0x4500 }};

    assert_ushort(H.v, ==, R.v);

    free(HX);

    return MUNIT_OK;
}
