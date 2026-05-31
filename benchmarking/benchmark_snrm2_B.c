//  Compares the classic running-scale snrm2 against the Anderson safe-scaling
//  snrm2_B on identical random vectors. Build (from repo root), e.g.:
//    cc -O2 -Iinclude -ISoftFloat/source/include benchmarking/benchmark_snrm2_B.c \
//       src/blas/level1/snrm2.c src/blas/level1/snrm2_B.c \
//       SoftFloat/build/<host>/softfloat.a -lm -o bench_snrm2_B
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "softblas.h"

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}

int main(void) {
    const int loops = 100;
    const int lengths[] = {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    const int nL = sizeof(lengths) / sizeof(lengths[0]);
    srand(0);  // fixed seed for reproducibility

    printf("%-9s %14s %14s %9s\n", "length", "snrm2 (s)", "snrm2_B (s)", "speedup");
    for (int li = 0; li < nL; li++) {
        int n = lengths[li];
        float32_t *SX = malloc(n * sizeof(float32_t));
        for (int i = 0; i < n; i++) {
            float r = 2.0f * rand() / RAND_MAX - 1.0f;   // values in [-1, 1]
            SX[i].v = *(uint32_t *)&r;
        }

        double b = now();
        for (int l = 0; l < loops; l++) { volatile float32_t d = snrm2(n, SX, 1, 'n'); (void)d; }
        double t_old = (now() - b) / loops;

        b = now();
        for (int l = 0; l < loops; l++) { volatile float32_t d = snrm2_B(n, SX, 1, 'n'); (void)d; }
        double t_new = (now() - b) / loops;

        printf("%-9d %14.6f %14.6f %8.2fx\n", n, t_old, t_new, t_old / t_new);
        free(SX);
    }
    return 0;
}
