#include <softblas.h>
#include "test.h"

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

//  Test SASUM, sum of absolute values.
//
//  float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX);
//
int main(int argc, char** argv) {
    union sing SX[] = {
        1.0f,
        2.0f,
        3.0f,
        4.0f,
        5.0f
    };
    uint64_t N = sizeof(SX) / sizeof(SX[0]);

    union sing assay = {.s = (float32_t) sasum(N, (float32_t*)SX, 1)};
    union sing canon = {.s = 15.0f};
    assert_int(assay.c, ==, canon.c);

    return 0;
}
