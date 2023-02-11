# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 26/02/23
 
import numpy as np

"""
Este programa calcula la integral de una funcion continua de una 
variable f(x) en el intervalo [a,b]
"""

def roomberg_00(a, b):
    return 0.5*(b-a)*(function(a)+function(b))

def romberg_n0(n, a, b):
    if (n==0):
        return roomberg_00(a, b)
    else:
        hn = (b-a)/(2**n)
        sum = 0
        for k in range (1, 2**(n-1) + 1):
            sum = sum + function(a + (2*k-1)*hn)
        Rn0 = 0.5*(romberg_n0(n-1, a, b)) + hn*sum
    return Rn0

def romberg_nm(n, m, a, b):
    # Definimos una tolarancia al error con cada iteración
    tol = 0.00001
    # Creamos una matríz donde almacenaremos las iteraciones de romberg
    matrix = np.zeros([n,m])
    matrix[0,0] = roomberg_00(a,b)
    
    for i in range(1, n):
        matrix[i,0] = romberg_n0(i, a, b)
        # Evaluamos el criterio de convergencia
        if (abs(matrix[i,0]-matrix[i-1,0]) <= tol ):
            print("La integral convergio")
            return matrix[i,0]

    for j in range(1, m):
        for i in range(j, n):
            aux = (4**(j))*matrix[i, j-1] - matrix[i-1, j-1]
            romberg_ij = (1.0/((4.0**(j))-1))*(aux)
            matrix[i,j] = romberg_ij
            # Evaluamos el criterio de convergencia
            if (abs(matrix[i,j]-matrix[i-1,j]) <= tol ):
                print("La integral convergio.")
                return matrix[i,j]

    # Si no converge nos quedamos con el valor más preciso que tenemos.
    return matrix[n,m]

"""
Funcion a integrar
"""
def function(x):
    g = 9.8; u = 1800; mo = 160000; q = 2500
    return u*np.log((mo)/(mo-q*x))-g*x

"""
Intervalo de integracion
"""
a = 0 ; b = 30
n = 10 ; m = 4

"""
Prueba del metodo
"""
integral_nm = romberg_nm(n, m, a, b)
print("La aproximacion de romberg(",n,",",m,") es: ", integral_nm)