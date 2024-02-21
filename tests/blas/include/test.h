#include <softfloat.h>
#include <softblas.h>

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

union half {
    float16_t s;
    uint16_t c;
};

union sing {
    float32_t s;
    uint32_t c;
};

union doub {
    float64_t s;
    uint64_t c;
};

union quad {
    float128_t* s;
    uint64_t* c;
};
