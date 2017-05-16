import numpy as np

def cholesky_reduce(A):
    pivot = A[0, 0]
    b = np.mat(A[1:, 0])
    B = A[1:, 1:]
    return B - (b.T * b) / pivot

def L(A):
    n = A.shape[0]
    if n == 1:
        return np.sqrt(A)
    b = np.mat(A[1:, 0])
    pivot = np.sqrt(A[0, 0])
    return np.bmat([
        [np.mat(pivot), np.zeros((1, n - 1))],
        [b.T / pivot, L(cholesky_reduce(A))]
    ])

def __main():
    A = np.array([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])
    print(L(A))

if __name__ == '__main__':
    __main()