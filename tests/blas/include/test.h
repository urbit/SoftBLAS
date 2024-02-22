#ifndef TEST_H
#define TEST_H

#include <softfloat.h>
#include <softblas.h>

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

union half {
    float16_t h;
    uint16_t u;
};

union sing {
    float32_t s;
    uint32_t u;
};

union doub {
    float64_t d;
    uint64_t u;
};

union quad {
    float128_t* q;
    uint64_t* u;
};

//  Test function prototypes
MunitResult test_sasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_sasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_saxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_saxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_saxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_saxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);

#endif // TEST_H
