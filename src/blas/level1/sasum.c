#include "softblas.h"

// typedef struct { uint32_t v; } float32_t;
// #define f32_abs(x) (float32_t)( (uint32_t)(x) & 0x7fffffff )

float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX) {
    float32_t stemp = { SB_REAL32_ZERO };
//     *(uint32_t*)(&stemp) = SB_REAL32_ZERO;
    
    uint64_t NincX = N * incX;
    for (uint64_t i = 0; i < NincX; i += incX) {
      // f32_abs(SX[i]);
      // ((uint32_t)SX[i].v) & 0x7fffffff;
      stemp = f32_add(stemp, f32_abs(SX[i]));
    }

    return nan_unify_s(stemp);
}
