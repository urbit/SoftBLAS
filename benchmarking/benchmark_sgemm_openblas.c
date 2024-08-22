#include <stdio.h>
#include <cblas.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main()
{
    const int num_loops_per_length = 1000;
    const int num_lengths = 9;
    const int lengths[] = {10, 20, 30, 40, 50, 75, 100, 200, 500};
    srand(0); // Fixed Seed for reproducibility
    for (int length_index = 0; length_index < num_lengths; length_index++)
    {
        int length = lengths[length_index];
        double total_time = 0.0;
        double total_time_squared = 0.0;
        for (int loop_num = 0; loop_num < num_loops_per_length; loop_num++)
        {
            float vector1[length * length];
            float vector2[length * length];
            float vector3[length * length];
            for (int vector_index = 0; vector_index < length; vector_index++)
            {
                vector1[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                vector2[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
                vector3[vector_index] = 2.0 * rand() / RAND_MAX - 1.0;
            }
            float* A = (float *)malloc(length * length * sizeof(float));
            memcpy(A, vector1, length * length * sizeof(float));
            float* B = (float *)malloc(length * length * sizeof(float));
            memcpy(B, vector2, length * length * sizeof(float));
            float* C = (float *)malloc(length * length * sizeof(float));
            memcpy(C, vector3, length * length * sizeof(float));
            float alpha = 2.0 * rand() / RAND_MAX - 1.0;
            float beta = 2.0 * rand() / RAND_MAX - 1.0;
            struct timespec begin, end;
            clock_gettime(CLOCK_MONOTONIC_RAW, &begin);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, length, length, length, alpha, A, length, B, length, beta, C, length);
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            double time_taken = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec - begin.tv_sec);
            total_time += time_taken;
            total_time_squared += time_taken * time_taken;
            free(A);
            free(B);
            free(C);
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