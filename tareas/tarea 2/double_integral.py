# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 14/02/23
 
import numpy as np


"""
Parametros de integracion
"""

# Definimos el integrando
def f(x, y):
    return np.log(x+2*y)

# definimos la region de integracion
xo = 1.4 ; xf = 2.0
yo = 1.0 ; yf = 1.5

# definimos el tamaño de subintervalo de integracion
dx = 0.0001 ; dy = 0.0001

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

def integracion_doble(xo, xf, yo, yf, dx, dy):
    Nx = int((xf-xo)/dx + 1)
    Ny = int((yf-yo)/dy + 1)
    suma = 0.0
    for i in range(1, Nx):
        suma = suma + trapecio_y(yo, yf, dy, xo + i*dx)
    aux = trapecio_y(yo,yf,dy,xo) + trapecio_y(yo,yf,dy,xf)
    integral = (0.5*dx)*(aux + 2*suma)
    return integral

"""
Mandamos a llamar al metodo
"""

integral = integracion_doble(xo, xf, yo, yf, dx, dy)
print("El valor de la integral es: ", integral)