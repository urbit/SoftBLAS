#include "test.h"

MunitResult test_hdot_0(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);

    float16_t HZ = hdot(5, (float16_t*)HX, 1, HY, 1);
    float16_t RZ = {*(uint16_t*)&(uint16_t){0x0 }};

    assert_ushort(HZ.v, ==, RZ.v);

    free(HX);
    free(HY);

    return MUNIT_OK;
}

MunitResult test_hdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x4500 , 0x4400 , 0x4200 , 0x4000 , 0x3c00 }, 5);

    float16_t HZ = hdot(5, (float16_t*)HX, 1, HY, 1);
    float16_t RZ = {*(uint16_t*)&(uint16_t){0x5060 }};

    assert_ushort(HZ.v, ==, RZ.v);

    free(HX);
    free(HY);

    return MUNIT_OK;
}

MunitResult test_hdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 }, 9);

    float16_t HZ = hdot(5, (float16_t*)HX, 1, HY, 2);
    float16_t RZ = {*(uint16_t*)&(uint16_t){0x4b80 }};

    assert_ushort(HZ.v, ==, RZ.v);

    free(HX);
    free(HY);

    return MUNIT_OK;
}

MunitResult test_hdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x4200 , 0x4400 , 0x4500 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x4500 , 0x4400 , 0x4200 , 0x4000 , 0x3c00 }, 5);

    float16_t HZ = hdot(5, (float16_t*)HX, 1, HY, -1);
    float16_t RZ = {*(uint16_t*)&(uint16_t){0x52e0 }};

    assert_ushort(HZ.v, ==, RZ.v);

    free(HX);
    free(HY);

    return MUNIT_OK;
}