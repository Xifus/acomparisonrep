import chaospy as cp
import numpy as np
from sobol import i4_sobol

#Add different sampling method

def create_sample_matrix(sample_matrix):
    N, k = sample_matrix.shape
    k = k/2
    A = sample_matrix[0:N, 0:k]
    B = sample_matrix[0:N, k:(2*k)]
    X = np.zeros(((2 + k) * N, k))
    X[0:N, 0:k] = A
    X[N:(2 * N), 0:k] = B
    for i in range(k):
        X[((2 + i) * N):((3 + i) * N), 0:k] = A
        X[((2 + i) * N):((3 + i) * N), i:(i + 1)] = B[0:N, i:(i + 1)]
    return X

def choose_sampling_method(N, k, rl):
    x = np.zeros((N, k))
    if rl == "R":
        x = np.random.random((N, k))
    elif rl == "S":
        x = cp.dist.sobol_lib.sobol(k, N).T
    elif rl == "L":
        x = cp.latin_hypercube(k, N).T
    elif rl == "sobol":
        seed = 123
        for i in range(N):
            r, seed = i4_sobol(k, seed)
            x[i, :] = r
    return x

def create_coefficient_a(k, a_type = "C"):
    a_array = np.zeros((k))
    if a_type == "A1-1":
        for i in range(2, k):
            a_array[i] = 6.52
    elif a_type == "A1-2":
        for i in range(0, k-2):
            a_array[i] = 6.52
    elif a_type == "A2":
        for i in range(5, k):
            a_array[i] = 6.52
    elif a_type == "B":
        for i in range(k):
            a_array[i] = 6.52
    else:
        a_array = a_array
    return a_array
        