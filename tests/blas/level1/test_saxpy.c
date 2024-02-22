#include "test.h"

//  Test SAXPY, A*X+Y.
//
//  void saxpy(uint64_t N, float32_t SA, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
//
MunitResult test_saxpy_0(const MunitParameter params[],
                         void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){0.0f}};
    float32_t SX[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){-2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){-4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };
    float32_t SY[] = {
        {*(uint32_t*)&(float){10.0f}},
        {*(uint32_t*)&(float){20.0f}},
        {*(uint32_t*)&(float){30.0f}},
        {*(uint32_t*)&(float){40.0f}},
        {*(uint32_t*)&(float){50.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    float32_t RY[] = {
        {*(uint32_t*)&(float){10.0f}},
        {*(uint32_t*)&(float){20.0f}},
        {*(uint32_t*)&(float){30.0f}},
        {*(uint32_t*)&(float){40.0f}},
        {*(uint32_t*)&(float){50.0f}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    return MUNIT_OK;
}

MunitResult test_saxpy_sum(const MunitParameter params[],
                           void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t SX[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){-2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){-4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };
    float32_t SY[] = {
        {*(uint32_t*)&(float){10.0f}},
        {*(uint32_t*)&(float){20.0f}},
        {*(uint32_t*)&(float){30.0f}},
        {*(uint32_t*)&(float){40.0f}},
        {*(uint32_t*)&(float){50.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    float32_t RY[] = {
        {*(uint32_t*)&(float){12.0f}},
        {*(uint32_t*)&(float){16.0f}},
        {*(uint32_t*)&(float){36.0f}},
        {*(uint32_t*)&(float){32.0f}},
        {*(uint32_t*)&(float){60.0f}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    return MUNIT_OK;
}

MunitResult test_saxpy_stride(const MunitParameter params[],
                           void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t SX[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };
    float32_t SY[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){1.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    float32_t RY[] = {
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){5.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){7.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){9.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){11.0f}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, 2);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    return MUNIT_OK;
}

MunitResult test_saxpy_neg_stride(const MunitParameter params[],
                           void *user_data) {
    const float32_t SA = {*(uint32_t*)&(float){2.0f}};
    float32_t SX[] = {
        {*(uint32_t*)&(float){1.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){4.0f}},
        {*(uint32_t*)&(float){5.0f}}
    };
    float32_t SY[] = {
        {*(uint32_t*)&(float){5.0f}},
        {*(uint32_t*)&(float){4.0f}},
        {*(uint32_t*)&(float){3.0f}},
        {*(uint32_t*)&(float){2.0f}},
        {*(uint32_t*)&(float){1.0f}}
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);
    float32_t RY[] = {
        {*(uint32_t*)&(float){15.0f}},
        {*(uint32_t*)&(float){12.0f}},
        {*(uint32_t*)&(float){9.0f}},
        {*(uint32_t*)&(float){6.0f}},
        {*(uint32_t*)&(float){3.0f}}
    };

    saxpy(N, (float32_t)SA, (float32_t*)SX, 1, (float32_t*)SY, -1);

    for (uint64_t i = 0; i < N; i++) {
        assert_int(SY[i].v, ==, RY[i].v);
    }

    return MUNIT_OK;
}
