#include "test.h"

MunitResult test_dswap_two(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);

    dswap(5, DX, 1, DY, 1);

    float64_t* RX = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);
    float64_t* RY = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DX[i].v, ==, RX[i].v);
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_dswap_stride(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 9);

    dswap(5, DX, 1, DY, 2);

    float64_t* RX = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);
    float64_t* RY = dvec((double[]){1.0, 0.0, -2.0, 0.0, 3.0, 0.0, -4.0, 0.0, 5.0}, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DX[i].v, ==, RX[i].v);
    }
    for (uint64_t i = 0; i < 9; i++) {
        assert_ullong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RX);
    free(RY);

    return MUNIT_OK;
}
