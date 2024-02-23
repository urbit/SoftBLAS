#include "test.h"

MunitResult test_haxpy_0(const MunitParameter params[], void *user_data) {
    const uint64_t N = 5;
    const float16_t HA = {*(uint16_t*)&(uint16_t){0.0}};
    float16_t* HX = hvec((uint16_t[]){0x3c00, 0xc000, 0x4000, 0xc000, 0x4800}, N);
    float16_t* HY = hvec((uint16_t[]){0x4120, 0x4160, 0x41a0, 0x41e0, 0x4220}, N);
    float16_t* RY = hvec((uint16_t[]){0x4120, 0x4160, 0x41a0, 0x41e0, 0x4220}, N);

    haxpy(N, HA, HX, 1, HY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_haxpy_sum(const MunitParameter params[],
                           void *user_data) {
    const uint64_t N = 5;
    const float16_t HA = {(uint16_t)0x4000};
    float16_t* HX = hvec((uint16_t[]){0x3c00, 0xc000, 0x4000, 0xc000, 0x4800}, N);
    float16_t* HY = hvec((uint16_t[]){0x4900, 0x4d00, 0x4f80, 0x5100, 0x5240}, N);
    float16_t* RY = hvec((uint16_t[]){0x4a00, 0x4c00, 0x5040, 0x5080, 0x5420}, N);

    haxpy(N, HA, HX, 1, HY, 1);
    
    for (uint64_t i = 0; i < N; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_haxpy_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float16_t HA = {(uint16_t)0x4000};
    float16_t* HX = hvec((uint16_t[]){0x3c00, 0xc000, 0x4000, 0xc400, 0x4800}, N);
    float16_t* HY = hvec((uint16_t[]){0x3c00, 0x4000, 0x4800, 0x4000, 0x4800, 0x4000, 0x4800, 0x4000, 0x4800}, 9);
    float16_t* RY = hvec((uint16_t[]){0x4200, 0x4000, 0x4400, 0x4000, 0x4a00, 0x4000, 0x0000, 0x4000, 0x4e00}, 9);

    haxpy(N, HA, HX, 1, HY, 2);

    for (uint64_t i = 0; i < 9; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_haxpy_neg_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float16_t HA = {(uint16_t)0x4000};
    float16_t* HX = hvec((uint16_t[]){0x3c00, 0x4000, 0x4200, 0x4400, 0x4500}, N);
    float16_t* HY = hvec((uint16_t[]){0x4500, 0x4400, 0x4200, 0x4000, 0x3c00}, N);
    float16_t* RY = hvec((uint16_t[]){0x4b80, 0x4a00, 0x4880, 0x4600, 0x4200}, N);

    haxpy(N, HA, HX, 1, HY, -1);

    for (uint64_t i = 0; i < N; i++) {
        assert_ushort(HY[i].v, ==, RY[i].v);
    }

    free(HX);
    free(HY);
    free(RY);

    return MUNIT_OK;
}
