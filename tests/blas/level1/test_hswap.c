#include "test.h"

MunitResult test_hswap_two(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);

    hswap(5, HX, 1, HY, 1);

    float16_t* RX = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);
    float16_t* RY = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HX[i].v, ==, RX[i].v);
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_hswap_stride(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);
    float16_t* HY = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 9);

    hswap(5, HX, 1, HY, 2);

    float16_t* RX = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);
    float16_t* RY = hvec((uint16_t[]){0x3c00 , 0x0 , 0xc000 , 0x0 , 0x4200 , 0x0 , 0xc400 , 0x0 , 0x4500 }, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HX[i].v, ==, RX[i].v);
    }
    for (uint64_t i = 0; i < 9; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}
