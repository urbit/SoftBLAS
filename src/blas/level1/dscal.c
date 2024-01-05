#include "softblas.h"

void dscal(uint64_t N,
           float64_t DA,
           float64_t *DX,
           uint64_t incX
          ) {
      if (incX == 1) {
            // code for increment equal to 1

            // clean-up loop
            uint64_t m = N % 5;
            if (m != 0) {
                  for (i = 0; i < m; i++) {
                        DX[i] = f64_mul(DA, DX[i]);
                  }
                  if (N < 5) return;
            }
            uint64_t mp1 = m + 1;
            for (uint64_t i = mp1; i < N; i += 5) {
                  DX[i] = f64_mul(DA, DX[i]);
                  DX[i + 1] = f64_mul(DA, DX[i + 1]);
                  DX[i + 2] = f64_mul(DA, DX[i + 2]);
                  DX[i + 3] = f64_mul(DA, DX[i + 3]);
                  DX[i + 4] = f64_mul(DA, DX[i + 4]);
            }
      } else {
            // code for increment not equal to 1
            uint64_t NincX = N * incX;
            for (uint64_t i = 0; i < NincX; i += incX) {
                  DX[i] = f64_mul(DA, DX[i]);
            }
      }
}
