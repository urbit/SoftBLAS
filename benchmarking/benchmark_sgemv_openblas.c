#include <stdio.h>
#include <cblas.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main()
{
    const int num_loops_per_length = 500;
    const int num_lengths = 7;
    const int lengths[] = {10, 20, 50, 100, 200, 500, 1000};
    srand(0); // Fixed Seed for reproducibility
    for (int length_index = 0; length_index < num_lengths; length_index++)
    {
        int length = lengths[length_index];
        double total_time = 0.0;
        double total_time_squared = 0.0;
        for (int loop_num = 0; loop_num < num_loops_per_length; loop_num++)
        {
            float vector1[length * length];
            float vector2[length];
            float vector3[length];
            for (int vector_index = 0; vector_index < length; vector_index++)
            {
                vector1[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                if (vector_index < length){
                    vector2[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                    vector3[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                }
            }
            float* A = (float *)malloc(length * length * sizeof(float));
            memcpy(A, vector1, length * length * sizeof(float));
            float* SX = (float *)malloc(length * sizeof(float));
            memcpy(SX, vector2, length * sizeof(float));
            float* SY = (float *)malloc(length * sizeof(float));
            memcpy(SY, vector3, length * sizeof(float));
            float alpha = 2.0 * rand() / RAND_MAX - 1.0;
            float beta = 2.0 * rand() / RAND_MAX - 1.0;
            struct timespec begin, end;
            clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
            cblas_sgemv(CblasRowMajor, CblasNoTrans, length, length, alpha, A, length, SX, 1, beta, SY, 1);
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            double time_taken = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec - begin.tv_sec);
            total_time += time_taken;
            total_time_squared += time_taken * time_taken;
            free(A);
            free(SX);
            free(SY);
        }
        printf(
            "Mean time taken for matrix of length %d: %.10f seconds +- %.10f seconds.\n",
            length, 
            total_time / num_loops_per_length, 
            sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)
        );
    }
    return 0;
}