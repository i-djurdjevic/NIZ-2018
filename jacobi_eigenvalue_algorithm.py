from numpy import array, identity, diagonal
from math import sqrt

#*****************************************************************************
#   Jakobijev algoritam za pronalazenje sopstvenih vrednosti matrica.
#
#   Ulazni parametri:
#     - n - velicina matrice
#     - a - kvadratna, realna i simetricna matrica dimenzije (n x n)
#     - it_max - maksimalan broj iteracija
#
#   Izlazni parametri:
#     - d - niz koji sadrzi sopstvene vrednosti matrice
#     - v - matrica koja sadrzi sopstvene vektore matrice
#     - it_num - ukupan broj iteracija
#     - rot_num - ukupan broj rotacija
#*****************************************************************************

def jacobi_eigenvalue_algorithm(a, n, tol = 1.0e-10):

    def max_off_diagonal(a, n): # Za pronalazenje najveceg vandijagonalnog elementa
        max_elem = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i,j]) >= max_elem:
                    max_elem = abs(a[i,j])
                    k = i
                    l = j
        return max_elem, k, l

    def rotate(a, n, p, k, l): # Rotiramo kako bismo anulirali a[k, l]
        diff = a[l,l] - a[k,k]

        if abs(a[k,l]) < abs(diff) * 0.1:
        	t = a[k,l] / diff
        else:
            theta = diff / (2.0 * a[k,l])
            t = 1.0 / (abs(theta) + sqrt(theta**2 + 1.0))
            
            if theta < 0.0:
            	t = -t

        c = 1.0 / sqrt(t**2 + 1.0);
        s = t*c
        tau = s / (1.0 + c)

        temp = a[k,l]
        a[k,l] = 0.0
        a[k,k] = a[k,k] - t * temp
        a[l,l] = a[l,l] + t * temp

        for i in range(k):      # i < k
            temp = a[i,k]
            a[i,k] = temp - s * (a[i,l] + tau * temp)
            a[i,l] = a[i,l] + s * (temp - tau * a[i,l])

        for i in range(k + 1, l):  # k < i < l
            temp = a[k,i]
            a[k,i] = temp - s * (a[i,l] + tau * a[k,i])
            a[i,l] = a[i,l] + s * (temp - tau * a[i,l])

        for i in range(l + 1, n):  # i > l
            temp = a[k,i]
            a[k,i] = temp - s * (a[l,i] + tau * temp)
            a[l,i] = a[l,i] + s * (temp - tau * a[l,i])

        for i in range(n):      # Azuriramo matricu transformacije
            temp = p[i,k]
            p[i,k] = temp - s * (p[i,l] + tau * p[i,k])
            p[i,l] = p[i,l] + s * (temp - tau * p[i,l])
        
    n = len(a)
    max_iter = 5 * (n**2)
    p = identity(n) * 1.0 

    for i in range(max_iter): 
        max_elem, k, l = max_off_diagonal(a, n)
        if max_elem < tol: 
        	return diagonal(a), p
        rotate(a, n, p, k, l)
    print ('Metoda nije konvergirala')

print(jacobi_eigenvalue_algorithm(array([[4, -30, 60, -35], [-30, 300, -675, 420], [60, -675, 1620, -1050], [-35, 420, -1050, 700]]), 4))
