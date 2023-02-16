! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 09/02/23

! En este programa resolvemos integrales dobles con la regla
! del trapecio.

        ! definimos la funcion que limita la region
        ! de integracion superior 
        real function limit_sup(x)
            implicit none
            real :: x
            limit_sup = 3.0 + exp(x/5.0)
            return
        end function

        ! definimos la funcion que limita la region
        ! de integracion inferior 
        real function limit_inf(x)
            implicit none
            real :: x
            limit_inf = log(x)
            return
        end function

        ! definimos el integrandpo 
        real function f(x,y)
            implicit none
            real :: x,y
            f = sin(x+y)
            return
        end function
 
        ! programa principal
        program main
 
            implicit none
            real :: integral, dx, dy, xo, xf

            ! definimos la región de integración
            xo = 1.0 ; xf = 3.0


            ! definimos el tamaño de subintervalo de integración
            dx = 0.001 ; dy = 0.001

            call integracion_doble(xo, xf, dx, dy, integral)

            write(*,*) "El valor de la integral es: ", integral
 
        end program main

        ! subrutina que cálcula la integral doble con regla del trapecio
        subroutine integracion_doble(xo, xf, dx, dy, integral)
            implicit none
            real :: xo, xf, yo, yf, dx, dy, integral, suma
            real :: limit_inf, limit_sup
            real :: trapecio_y, aux
            integer :: Nx, i
            Nx = (xf-xo)/dx + 1
            suma = 0.0
            do i = 1, Nx-1
                yo =  limit_inf(xo + i*dx)
                yf = limit_sup(xo + i*dx)
                suma = suma + trapecio_y(yo, yf, dy, xo + i*dx)
            end do
            aux = trapecio_y(limit_inf(xo),limit_sup(xo),dy,xo) 
            aux = aux + trapecio_y(limit_inf(xo),limit_sup(xo),dy,xf)
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