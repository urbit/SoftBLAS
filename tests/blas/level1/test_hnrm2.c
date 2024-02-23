#include "test.h"

MunitResult test_hnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x0 , 0x0 }, 2);

    float16_t H = (float16_t) hnrm2(2, (float16_t*)HX, 1);
    float16_t R = (float16_t){0x0 };

    assert_ushort(H.v, ==, R.v);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_hnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);

    float16_t H = (float16_t) hnrm2(5, (float16_t*)HX, 1);
    float16_t R = {*(uint16_t*)&(uint16_t){0x476b}};  // sqrt(55.0)

    assert_ushort(H.v, ==, R.v);
    
    free(HX);

    return MUNIT_OK;
}

MunitResult test_hnrm2_stride(const MunitParameter params[],
                              void *user_data) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 }, 9);

    float16_t H = (float16_t) hnrm2(5, (float16_t*)HX, 2);
    float16_t R = {*(uint16_t*)&(uint16_t){0x4079}};  // sqrt(5.0)

    assert_ushort(H.v, ==, R.v);

    free(HX);

    return MUNIT_OK;
}
