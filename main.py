import numpy as np
from termcolor import colored
# Name 1: Tehila Levi ID: 322257213
# Name 2: Shani Sapir ID: 323060590


def find_det(a):
    """"" Returns The Determinant Of A Matrix"""
    s = np.array(a)
    return np.linalg.det(s)


def copy_m(a, size):
    new_mat = []
    for i in range(size):
        row_list = []
        for j in range(size):
            row_list.append(a[i][j])
        new_mat.append(row_list)
    return new_mat


def is_dominant(a, n):
    mat = copy_m(a, n)
    for i in range(0, n):
        _sum = 0
        for j in range(0, n):
            _sum = _sum + abs(mat[i][j])
        _sum = _sum - abs(mat[i][i])
        if abs(a[i][i]) < _sum:
            return False
    return True


def initialize_d(a, n):
    new_mat = copy_m(a, n)
    for i in range(n):
        val = abs(new_mat[i][i])
        for j in range(n):
            (new_mat[i][j]) = 0
        (new_mat[i][i]) = val
    return new_mat


def calc_m1(a, b):
    """"" Returns The Result Matrix Of 2 Matrix Calculation"""

    matrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):
                matrix[i][j] += a[i][k] * b[k][j]
    return matrix


def calc_m2(a, b):
    """"" Calculates 3x1 matrix"""
    new_mat = [[0, ], [0, ], [0, ]]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(a)):
                new_mat[i][j] += (a[i][k] * b[k][j])

    return new_mat


def inverse_matrix(a):
    """"" Calculating An Inverse Of 4X4 Matrix According To Gauss Elimination Method Using Elementary Matrices.
    """
    if find_det(a) != 0:
        i = np.identity(len(a))
        i_copy = i
        a_copy = a
        indices = list(range(len(a)))  # indices = [0, 1, 2, 3]
        for i in range(len(a)):  # any element in range 0 - 3
            div = 1.0/a_copy[i][i]
            for j in range(len(a)):
                a_copy[i][j] *= div
                i_copy[i][j] *= div
            for s in indices[0:i]+indices[i+1:]:
                # from index 0 until i (not including i) and from index i+1 until the end
                val = a_copy[s][i]
                for x in range(len(a)):
                    a_copy[s][x] -= val*a_copy[i][x]
                    i_copy[s][x] -= val*i_copy[i][x]
    return i_copy


def lu_de_comp(a, n):
    l_ = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u_ = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(len(a)):
        for j in range(len(a)):
            if i != j:
                if i > j:
                    l_[i][j] = a[i][j]
                else:
                    u_[i][j] = a[i][j]

    return l_, u_


def negative_m(a):
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] = -a[i][j]
    return a


def sum_mat(a, b):
    for i in range(len(a)):
        for j in range(len(a[0])):
            a[i][j] += b[i][j]
    return a


def lu_solution(mat, n):
    l, u = lu_de_comp(mat, n)
    d = initialize_d(mat, n)
    d1 = inverse_matrix(d)
    return l, u, d, d1


def sub_mat(a, b):
    sub = 0
    for i in range(len(a)):
        for j in range(1):
            sub += a[i][j] - b[i][j]
    return sub


def gauss_seidel_method(a, n, b):
    """"" Gets  NxN Matrix and returns Calculation Approximation According To Gauss Seidel Method"""
    l, u = lu_de_comp(a, n)
    d = [[4, 0, 0], [0, 10, 0], [0, 0, 5]]
    # diagonal
    xr = [[0], [0], [0]]
    xr1 = [[0], [0], [0]]
    sum_ld = sum_mat(l, d)  # (L+D)
    h = inverse_matrix(sum_ld)  # (L+D)^-1
    m = negative_m(h)  # -(L+D)^-1
    g = calc_m1(m, u)  # -(L+D)^-1 * U
    xr1 = sum_mat(calc_m2(g, xr), calc_m1(h, b))  # According to Formula
    mat = []
    xr2 = []
    while abs(sub_mat(xr1, xr) > 0.001):
        xr = xr1
        xr1 = sum_mat(calc_m2(g, xr), calc_m2(h, b))
        for i in xr1:
            z = [np.absolute(i) for i in xr1]
        for x, j in enumerate(z):
            print(j)
            if x % 3 == 2:
                print("\t")
    return


def yakobi_method(mat, n, b):
    xr = [[0], [0], [0]]
    xr1 = [[0], [0], [0]]
    l, u, d, d1 = lu_solution(mat, n)
    h = d1
    g = negative_m(calc_m1(d1, sum_mat(l, u)))
    xr1 = sum_mat(calc_m2(g, xr), calc_m2(h, b))
    while abs(sub_mat(xr1, xr)) > 0.001:
        xr = xr1
        xr1 = sum_mat(calc_m2(g, xr), calc_m2(h, b))
        for i in xr1:
            z = [np.absolute(i) for i in xr1]
        for x, j in enumerate(z):
            print(j)
            if x % 3 == 2:
                print("\t")
    return


def main():
    """"" Activating All Functions together"""

    a = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
    b = [[2], [6], [5]]
    n = len(a)
    d = initialize_d(a, n)
    if not is_dominant(a, n):
        print("No Dominant diagonal Was Found.")
    else:
        print(colored("Yakobi_Method:", "red", attrs=['bold']))
        yakobi_method(a, n, b)
        print(colored("Gauss_Seidel_Method:", "red", attrs=['bold']))
        print(gauss_seidel_method(a, n, b))


main()
