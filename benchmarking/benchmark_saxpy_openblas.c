#include <stdio.h>
#include <cblas.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main()
{
    const int num_loops_per_length = 100;
    const int num_lengths = 7;
    const int lengths[] = {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    srand(0); // Fixed Seed for reproducibility
    for (int length_index = 0; length_index < num_lengths; length_index++)
    {
        int length = lengths[length_index];
        double total_time = 0.0;
        double total_time_squared = 0.0;
        for (int loop_num = 0; loop_num < num_loops_per_length; loop_num++)
        {
            float vector1[length];
            float vector2[length];
            for (int vector_index = 0; vector_index < length; vector_index++)
            {
                vector1[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                vector2[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
            }
            float SA = 2.0 * rand() / RAND_MAX - 1.0;
            float *SX = (float *)malloc(length * sizeof(float));  // Use a Heap Array so that it's a fair comparison
            float *SY = (float *)malloc(length * sizeof(float));
            memcpy(SX, vector1, length * sizeof(float));
            memcpy(SY, vector2, length * sizeof(float));
            struct timespec begin, end;
            clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
            cblas_saxpy(length, SA, SX, 1, SY, 1);
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            double time_taken = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec - begin.tv_sec);
            total_time += time_taken;
            total_time_squared += time_taken * time_taken;
            free(SX);
            free(SY);
        }
        printf(
            "Mean time taken for vector of length %d: %.10f seconds +- %.10f seconds.\n",
            length, 
            total_time / num_loops_per_length, 
            sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)
        );
    }
    return 0;
}