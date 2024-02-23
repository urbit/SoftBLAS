#include "test.h"

MunitResult test_dscal_0(const MunitParameter params[], void *user_data) {
    const float64_t DA = {*(uint64_t*)&(double){0.0}};
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);

    dscal(5, DA, DX, 1);

    float64_t* RX = dvec((double[]){0.0, -0.0, 0.0, -0.0, 0.0}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DX[i].v, ==, RX[i].v);
    }

    free(DX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_dscal_12345(const MunitParameter params[],
                           void *user_data) {
    const float64_t DA = {*(uint64_t*)&(double){10.5}};
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);

    dscal(5, DA, DX, 1);

    float64_t* RX = dvec((double[]){10.5, -21.0, 31.5, -42.0, 52.5}, 5);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DX[i].v, ==, RX[i].v);
    }

    free(DX);
    free(RX);

    return MUNIT_OK;
}

MunitResult test_dscal_stride(const MunitParameter params[],
                              void *user_data) {
    const float64_t DA = {*(uint64_t*)&(double){2.5}};
    float64_t* DX = dvec((double[]){1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0}, 9);

    dscal(5, DA, DX, 2);

    float64_t* RX = dvec((double[]){2.5, 2.0, 2.5, 2.0, 2.5, 2.0, 2.5, 2.0, 2.5}, 9);

    for (uint64_t i = 0; i < 5; i++) {
        assert_ullong(DX[i].v, ==, RX[i].v);
    }

    free(DX);
    free(RX);

    return MUNIT_OK;
}
