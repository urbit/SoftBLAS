#include "test.h"

MunitResult test_hcopy_all(const MunitParameter params[],
                           void* user_data_or_fixture) {
    const uint64_t N = 5;
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0xc000 , 0x4200 , 0xc400 , 0x4500 }, N);
    float16_t* HY = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, N);

    hcopy(N, (float16_t*)HX, 1, HY, 1);

    float16_t RY[] = {
        {*(uint16_t*)&(uint16_t){0x3c00 }},
        {*(uint16_t*)&(uint16_t){0xc000 }},
        {*(uint16_t*)&(uint16_t){0x4200 }},
        {*(uint16_t*)&(uint16_t){0xc400 }},
        {*(uint16_t*)&(uint16_t){0x4500 }}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);

    return MUNIT_OK;
}

MunitResult test_hcopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float16_t* HX = hvec((uint16_t[]){0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 , 0x4000 , 0x3c00 }, 9);
    float16_t* HY = hvec((uint16_t[]){0x0 , 0x0 , 0x0 , 0x0 , 0x0 }, 5);

    hcopy(5, (float16_t*)HX, 2, HY, 1);

    float16_t RY[] = {
        {*(uint16_t*)&(uint16_t){0x3c00 }},
        {*(uint16_t*)&(uint16_t){0x3c00 }},
        {*(uint16_t*)&(uint16_t){0x3c00 }},
        {*(uint16_t*)&(uint16_t){0x3c00 }},
        {*(uint16_t*)&(uint16_t){0x3c00 }}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);

    return MUNIT_OK;
}
