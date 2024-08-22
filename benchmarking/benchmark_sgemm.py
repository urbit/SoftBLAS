from math import sqrt
from random import uniform, seed
from time import monotonic_ns

from numpy import ndarray, array, float32, dot

num_loops_per_length: int = 200
lengths: list[int] = [10, 20, 30, 40, 50, 75, 100, 200, 500]
seed(0)
for length in lengths:
    total_time: float = 0.0
    total_time_squared: float = 0.0
    for _ in range(num_loops_per_length):
        matrix1: ndarray = array([[uniform(-1, 1) for _ in range(length)] for _ in range(length)], dtype=float32)
        matrix2: ndarray = array([[uniform(-1, 1) for _ in range(length)] for _ in range(length)], dtype=float32)
        matrix3: ndarray = array([[uniform(-1, 1) for _ in range(length)] for _ in range(length)], dtype=float32)
        alpha: float = uniform(-1, 1)
        beta: float = uniform(-1, 1)
        start: int = monotonic_ns()
        matrix3 = alpha * matrix1 @ matrix2 + beta * matrix3
        end: int = monotonic_ns()
        total_time += (end - start) / 10 ** 9
        total_time_squared += (end - start) ** 2 / 10 ** 18
    print(
            f"Mean time taken for matrix of length {length}: {total_time / num_loops_per_length} seconds "
            f"+- {sqrt(total_time_squared / num_loops_per_length - (total_time / num_loops_per_length) * (total_time / num_loops_per_length)) / sqrt(num_loops_per_length)} seconds."
    )