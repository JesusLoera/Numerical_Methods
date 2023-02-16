# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 13/02/23
 
import numpy as np

"""
Parametros de integracion
"""

# Definimos el integrando
def f(x, y):
    return np.sin(x+y)

# Definimos el límite superior de integracion
def limit_sup(x):
    return 3 + np.exp(x/5.0)

# Definimos el límite superior de integracion
def limit_inf(x):
    return np.log(x)

# definimos la region de integracion
xo = 1.0 ; xf = 3.0

# definimos el tamaño de subintervalo de integracion
dx = 0.001 ; dy = 0.001

"""
Funciones que definen el metodo de integracion
"""

def trapecio_y(yo, yf, dy, xcte):
    N = int((yf-yo)/dy + 1)
    suma = 0.0
    for i in range(1, N):
        suma = suma + f(xcte, yo+i*dy)
    trapecio_y = (0.5*dy)*(f(xcte,yo) + f(xcte,yf) + 2*suma)
    return trapecio_y

def integracion_doble(xo, xf, dx, dy):
    Nx = int((xf-xo)/dx + 1)
    suma = 0.0
    for i in range(1, Nx):
        yo = limit_inf(xo)
        yf = limit_sup(xf)
        suma = suma + trapecio_y(yo, yf, dy, xo + i*dx)
    yo = limit_inf(xo)
    yf = limit_sup(xf)
    aux = trapecio_y(yo,yf,dy,xo) + trapecio_y(yo,yf,dy,xf)
    integral = (0.5*dx)*(aux + 2*suma)
    return integral

"""
Mandamos a llamar al metodo
"""

integral = integracion_doble(xo, xf, dx, dy)
print("El valor de la integral es: ", integral)