from math import sqrt
from random import uniform, seed
from time import monotonic_ns

from numpy import ndarray, array, float32, fabs, argmax

num_loops_per_length: int = 100
lengths: list[int] = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
seed(0)
for length in lengths:
    total_time: float = 0.0
    total_time_squared: float = 0.0
    for _ in range(num_loops_per_length):
        vector: ndarray = array([uniform(-1, 1) for _ in range(length)], dtype=float32)
        start: int = monotonic_ns()
        dummy: int = argmax(fabs(vector))
        end: int = monotonic_ns()
        total_time += (end - start) / 10 ** 9
        total_time_squared += (end - start) ** 2 / 10 ** 18
    print(
            f"Mean time taken for vector of length {length}: {total_time / num_loops_per_length} seconds "
            f"+- {sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)} seconds."
    )