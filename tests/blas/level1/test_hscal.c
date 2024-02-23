#include "test.h"

MunitResult test_hscal_0(const MunitParameter params[], void *user_data) {
    const float16_t HA = {*(uint16_t*)&(uint16_t){0x0 }};
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);

    hscal(5, HA, HX, 1);

    float16_t* RX = hvec((uint16_t[]){0x0 , 0x8000 , 0x0 , 0x8000 , 0x0 }, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HX[i].v, ==, RX[i].v);
    }

    free(HX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_hscal_12345(const MunitParameter params[],
                           void *user_data) {
    const float16_t HA = {*(uint16_t*)&(uint16_t){0x4940}};
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, 5);

    hscal(5, HA, HX, 1);

    float16_t* RX = hvec((uint16_t[]){0x4940, 0xcd40 , 0x4fe0, 0xd140 , 0x5290}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HX[i].v, ==, RX[i].v);
    }

    free(HX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_hscal_stride(const MunitParameter params[],
                              void *user_data) {
    const float16_t HA = {*(uint16_t*)&(uint16_t){0x4100}};
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 }, 9);

    hscal(5, HA, HX, 2);

    float16_t* RX = hvec((uint16_t[]){0x4100, 0x4000 , 0x4100, 0x4000 , 0x4100, 0x4000 , 0x4100, 0x4000 , 0x4100}, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HX[i].v, ==, RX[i].v);
    }

    free(HX);
    free(RX);

    return MUNIT_OK;
}