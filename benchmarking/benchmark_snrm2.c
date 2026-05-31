//  Times the classic running-scale snrm2 over a range of vector lengths.
//  Build (from repo root), e.g.:
//    cc -O2 -Iinclude -ISoftFloat/source/include benchmarking/benchmark_snrm2.c \
//       src/blas/level1/snrm2.c SoftFloat/build/<host>/softfloat.a -lm -o bench_snrm2
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "softblas.h"

int main()
{
    const int num_loops_per_length = 100;
    const int num_lengths = 7;
    const int lengths[] = {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    const int max_length = 1000000;  // guard: largest vector we ever allocate
    srand(0); // Fixed Seed for reproducibility

    //  Heap-allocate the working buffer once, up front, instead of using a
    //  stack VLA (`float32_t vector[length]`) per loop iteration. A VLA does
    //  not fail gracefully like malloc: an oversized length walks off the
    //  stack into mapped pages and thrashes RAM rather than returning NULL.
    float32_t *SX = malloc((size_t)max_length * sizeof(float32_t));
    if (SX == NULL) { perror("malloc"); return 1; }

    for (int length_index = 0; length_index < num_lengths; length_index++)
    {
        int length = lengths[length_index];
        if (length < 1 || length > max_length)
        {
            fprintf(stderr, "skipping out-of-range length %d (max %d)\n", length, max_length);
            continue;
        }
        double total_time = 0.0;
        double total_time_squared = 0.0;
        for (int loop_num = 0; loop_num < num_loops_per_length; loop_num++)
        {
            for (int vector_index = 0; vector_index < length; vector_index++)
            {
                float r = 2.0f * rand() / RAND_MAX - 1.0f;  // values in [-1, 1]
                SX[vector_index].v = *(uint32_t *)&r;
            }
            struct timespec begin, end;
            clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
            volatile float32_t dummy = snrm2(length, SX, 1, 'n');
            (void)dummy;
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            double time_taken = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec - begin.tv_sec);
            total_time += time_taken;
            total_time_squared += time_taken * time_taken;
        }
        printf(
            "Mean time taken for vector of length %d: %f seconds +- %f seconds.\n",
            length,
            total_time / num_loops_per_length,
            sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)
        );
    }
    free(SX);
    return 0;
}
