from numpy import array, identity, diagonal, zeros
from math import sqrt

# *****************************************************************************
#   Jakobijev algoritam za pronalazenje sopstvenih vrednosti matrica.
#
#   Ulazni parametri:
#     - n - velicina matrice
#     - a - kvadratna, realna i simetricna matrica dimenzije (n x n)
#     - max_iter - maksimalan broj iteracija
#
#   Izlazni parametri:
#     - diagonal(a) - niz koji sadrzi sopstvene vrednosti matrice
#     - p - matrica koja sadrzi sopstvene vektore matrice
#     - it_num - ukupan broj iteracija
# *****************************************************************************

def print_results(eigenvalues, eigenvectors, num_of_it):
    print('Eigenvalues:', eigenvalues)
    print('Eigenvectors:', eigenvectors)
    print('Total number of rotations:', num_of_it)
    print('\n')

def jacobi_eigenvalue_algorithm(a, n, max_iter):

    def treshold_convergence(a, n, it_num):
        off_sum = 0.0

        for j in range(0, n):
            for i in range(0, j):
                off_sum = off_sum + a[i, j] ** 2

        off_sum = sqrt(off_sum) / float(4 * n)

        if (off_sum == 0.0):
            return off_sum

        for p in range(0, n):
            for q in range(p + 1, n):
                gapq = 10.0 * abs(a[p, q])
                termp = gapq + abs(a[p, p])
                termq = gapq + abs(a[q, q])


        if 4 < it_num and termp == abs(a[p, p]) and termq == abs(a[q, q]):
            a[p, q] = 0.0

        return off_sum


    def max_off_diagonal(a, n):  # Za pronalazenje najveceg vandijagonalnog elementa
        max_elem = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= max_elem:
                    max_elem = abs(a[i, j])
                    k = i
                    l = j
        return max_elem, k, l


    def rotate(a, n, p, k, l):  # Rotiramo kako bismo anulirali a[k, l]
        diff = a[l, l] - a[k, k]

        if abs(a[k, l]) < abs(diff) * 0.1:
            t = a[k, l] / diff
        else:
            theta = diff / (2.0 * a[k, l])
            t = 1.0 / (abs(theta) + sqrt(theta ** 2 + 1.0))

            if theta < 0.0:
                t = -t

        c = 1.0 / sqrt(t ** 2 + 1.0);
        s = t * c
        tau = s / (1.0 + c)

        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp

        for i in range(k):  # i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])

        for i in range(k + 1, l):  # k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])

        for i in range(l + 1, n):  # i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])

        for i in range(n):  # Azuriramo matricu transformacije
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])


    n = len(a)
    p = identity(n) * 1.0
    for it_num in range(max_iter):
        max_elem, k, l = max_off_diagonal(a, n)
        off_sum = treshold_convergence(a, n, it_num)

        if max_elem < 0.000001 or off_sum == 0.0:
            return diagonal(a), p, it_num

        rotate(a, n, p, k, l)

    print('Metoda nije konvergirala')

# Tests:

eigenvalues, eigenvectors, num_of_it = jacobi_eigenvalue_algorithm(
    array([\
        [4, -30, 60, -35], \
        [-30, 300, -675, 420], \
        [60, -675, 1620, -1050], \
        [-35, 420, -1050, 700]]), 4, 100)

print_results(eigenvalues, eigenvectors, num_of_it)

eigenvalues, eigenvectors, num_of_it = jacobi_eigenvalue_algorithm(array([ \
    [ 4.0, 0.0, 0.0, 0.0 ], \
    [ 0.0, 1.0, 0.0, 0.0 ], \
    [ 0.0, 0.0, 3.0, 0.0 ], \
    [ 0.0, 0.0, 0.0, 2.0 ] ] ), 4, 100)

print_results(eigenvalues, eigenvectors, num_of_it)

n = 5
a = zeros([n, n])
for j in range ( 0, n ):
    for i in range ( 0, n ):
        if ( i == j ):
            a[i,j] = -2.0
        elif ( i == j + 1 or i == j - 1 ):
            a[i,j] = 1.0


eigenvalues, eigenvectors, num_of_it = jacobi_eigenvalue_algorithm(a, n, 100)

print_results(eigenvalues, eigenvectors, num_of_it)

