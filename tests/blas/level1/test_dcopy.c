#include "test.h"

MunitResult test_dcopy_all(const MunitParameter params[],
                           void* user_data_or_fixture) {
    const uint64_t N = 5;
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, N);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, N);

    dcopy(N, (float64_t*)DX, 1, DY, 1);

    float64_t RY[] = {
        {*(uint64_t*)&(double){1.0}},
        {*(uint64_t*)&(double){-2.0}},
        {*(uint64_t*)&(double){3.0}},
        {*(uint64_t*)&(double){-4.0}},
        {*(uint64_t*)&(double){5.0}}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);

    return MUNIT_OK;
}

MunitResult test_dcopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0}, 9);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);

    dcopy(5, (float64_t*)DX, 2, DY, 1);

    float64_t RY[] = {
        {*(uint64_t*)&(double){1.0}},
        {*(uint64_t*)&(double){1.0}},
        {*(uint64_t*)&(double){1.0}},
        {*(uint64_t*)&(double){1.0}},
        {*(uint64_t*)&(double){1.0}}
    };

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);

    return MUNIT_OK;
}
