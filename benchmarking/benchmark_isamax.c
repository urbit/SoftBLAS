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
    srand(0); // Fixed Seed for reproducibility
    for (int length_index = 0; length_index < num_lengths; length_index++)
    {
        int length = lengths[length_index];
        double total_time = 0.0;
        double total_time_squared = 0.0;
        for (int loop_num = 0; loop_num < num_loops_per_length; loop_num++)
        {
            float vector[length];
            for (int vector_index = 0; vector_index < length; vector_index++)
            {
                vector[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
            }
            float32_t *SX = svec(vector, length);
            struct timespec begin, end;
            clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
            volatile uint64_t dummy = isamax(length, SX, 1);
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            double time_taken = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec - begin.tv_sec);
            total_time += time_taken;
            total_time_squared += time_taken * time_taken;
            free(SX);
        }
        printf(
            "Mean time taken for vector of length %d: %f seconds +- %f seconds.\n",
            length, 
            total_time / num_loops_per_length, 
            sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)
        );
    }
    return 0;
}