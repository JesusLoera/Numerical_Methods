! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 09/02/23

! En este programa resolvemos integrales dobles con la regla
! del trapecio.

        ! definimos el integrandpo 
        real function f(x,y)
            implicit none
            real :: x,y
            f = log(x+2*y)
            return
        end function

        ! programa principal
        program main
 
            implicit none
            real :: integral, dx, dy, xo, xf, yo, yf

            ! definimos la región de integración
            xo = 1.4 ; xf = 2.0
            yo = 1.0 ; yf = 1.5

            ! definimos el tamaño de subintervalo de integración
            dx = 0.0001 ; dy = 0.0001

            call integracion_doble(xo, xf, yo, yf, dx, dy, integral)

            write(*,*) "El valor de la integral es: ", integral
 
        end program main

        ! subrutina que cálcula la integral doble con regla del trapecio
        subroutine integracion_doble(xo, xf, yo, yf, dx, dy, integral)
            implicit none
            real :: xo, xf, yo, yf, dx, dy, integral, suma
            real :: f
            real :: trapecio_y, aux
            integer :: Nx, Ny, i
            Nx = (xf-xo)/dx + 1
            Ny = (yf-yo)/dy + 1
            suma = 0.0
            do i = 1, Nx-1
                suma = suma + trapecio_y(yo, yf, dy, xo + i*dx)
            end do
            aux = trapecio_y(yo,yf,dy,xo) + trapecio_y(yo,yf,dy,xf)
            integral = (0.5*dx)*(aux + 2*suma)
        end subroutine

        ! regla del trapecio para integrar funciones f(cte,y)
        real function trapecio_y(yo, yf, dy, xcte)
            implicit none
            real :: yo, yf, dy, suma
            ! función del integrando
            real :: f, xcte
            integer :: N, i

            N = (yf-yo)/dy + 1
            suma = 0.0
            do i = 1, N-1
                suma = suma + f(xcte, yo+i*dy)
            end do
            trapecio_y = (0.5*dy)*(f(xcte,yo) + f(xcte,yf) + 2*suma)

            return
        end function