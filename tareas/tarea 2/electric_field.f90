! Programa elaborado por Jesus Eduardo Loera Casas
! Fecha 09/02/23

! En este programa resolvemos integrales dobles con la regla
! del trapecio.

        ! programa principal
        program main

            implicit none
            real :: dx, dy, xo, xf
            real :: Ex, Ey, Ez

            ! definimos la región de integración
            xo = -2.0 ; xf = +2.0

            ! definimos el tamaño de subintervalo de integración
            dx = 0.001 ; dy = 0.001

            call integracion_doble_fx(xo, xf, dx, dy, Ex)
            write(*,*) "La componente Ex es: ", Ex

            call integracion_doble_fy(xo, xf, dx, dy, Ey)
            write(*,*) "La componente Ey es: ", Ey

            call integracion_doble_fz(xo, xf, dx, dy, Ez)
            write(*,*) "La componente Ez es: ", Ez

        end program main

        ! definimos el integrando Ez
        real function fz(x,y)
            implicit none
            real :: x,y
            real :: pi, sigma, epsilon, zo
            pi = 3.14159265359
            epsilon = 8.85e-12
            sigma = 5.0/(4.0*pi)
            zo =  5
            fz = (sigma*zo)/( 4*pi*epsilon*(x**2 + y**2 + zo**2)**(1.5) )
            return
        end function

        ! definimos el integrando para Ex
        real function fx(x,y)
            implicit none
            real :: x,y
            real :: pi, sigma, epsilon, zo
            pi = 3.14159265359
            epsilon = 8.85e-12
            sigma = 5.0/(4.0*pi)
            zo =  5
            fx = -(sigma*x)/( 4*pi*epsilon*(x**2 + y**2 + zo**2)**(1.5) )
            return
        end function

        ! definimos el integrando para Ey
        real function fy(x,y)
            implicit none
            real :: x,y
            real :: pi, sigma, epsilon, zo
            pi = 3.14159265359
            epsilon = 8.85e-12
            sigma = 5.0/(4.0*pi)
            zo =  5
            fy = -(sigma*y)/( 4*pi*epsilon*(x**2 + y**2 + zo**2)**(1.5) )
            return
        end function

        ! definimos la funcion que limita la region
        ! de integracion superior 
        real function limit_sup(x)
            implicit none
            real :: x
            limit_sup = sqrt(4.0-x**2)
            return
        end function

        ! definimos la funcion que limita la region
        ! de integracion inferior 
        real function limit_inf(x)
            implicit none
            real :: x
            limit_inf = -sqrt(4.0-x**2)
            return
        end function

        ! subrutina que cálcula la integral doble con regla del trapecio
        subroutine integracion_doble_fx(xo, xf, dx, dy, integral)
            implicit none
            real :: xo, xf, yo, yf, dx, dy, integral, suma
            real :: limit_inf, limit_sup
            real :: trapeciofx_y, aux
            integer :: Nx, i
            Nx = (xf-xo)/dx + 1
            suma = 0.0
            do i = 1, Nx-1
                yo =  limit_inf(xo + i*dx)
                yf = limit_sup(xo + i*dx)
                suma = suma + trapeciofx_y(yo, yf, dy, xo + i*dx)
            end do
            aux = trapeciofx_y(limit_inf(xo),limit_sup(xo),dy,xo) 
            aux = aux + trapeciofx_y(limit_inf(xo),limit_sup(xo),dy,xf)
            integral = (0.5*dx)*(aux + 2*suma)
        end subroutine

        ! subrutina que cálcula la integral doble con regla del trapecio
        subroutine integracion_doble_fy(xo, xf, dx, dy, integral)
            implicit none
            real :: xo, xf, yo, yf, dx, dy, integral, suma
            real ::  limit_inf, limit_sup
            real :: trapeciofy_y, aux
            integer :: Nx, i
            Nx = (xf-xo)/dx + 1
            suma = 0.0
            do i = 1, Nx-1
                yo =  limit_inf(xo + i*dx)
                yf = limit_sup(xo + i*dx)
                suma = suma + trapeciofy_y(yo, yf, dy, xo + i*dx)
            end do
            aux = trapeciofy_y(limit_inf(xo),limit_sup(xo),dy,xo) 
            aux = aux + trapeciofy_y(limit_inf(xo),limit_sup(xo),dy,xf)
            integral = (0.5*dx)*(aux + 2*suma)
        end subroutine

        ! subrutina que cálcula la integral doble con regla del trapecio
        subroutine integracion_doble_fz(xo, xf, dx, dy, integral)
            implicit none
            real :: xo, xf, yo, yf, dx, dy, integral, suma
            real ::  limit_inf, limit_sup
            real :: trapeciofz_y, aux
            integer :: Nx, i
            Nx = (xf-xo)/dx + 1
            suma = 0.0
            do i = 1, Nx-1
                yo =  limit_inf(xo + i*dx)
                yf = limit_sup(xo + i*dx)
                suma = suma + trapeciofz_y(yo, yf, dy, xo + i*dx)
            end do
            aux = trapeciofz_y(limit_inf(xo),limit_sup(xo),dy,xo) 
            aux = aux + trapeciofz_y(limit_inf(xo),limit_sup(xo),dy,xf)
            integral = (0.5*dx)*(aux + 2*suma)
        end subroutine

        ! regla del trapecio para integrar funciones fx(cte,y)
        real function trapeciofx_y(yo, yf, dy, xcte)
            implicit none
            real :: yo, yf, dy, suma
            ! función del integrando
            real :: fx, xcte
            integer :: N, i

            N = (yf-yo)/dy + 1
            suma = 0.0
            do i = 1, N-1
                suma = suma + fx(xcte, yo+i*dy)
            end do
            trapeciofx_y = (0.5*dy)*(fx(xcte,yo) + fx(xcte,yf) + 2*suma)
            return
        end function

        ! regla del trapecio para integrar funciones fy(cte,y)
        real function trapeciofy_y(yo, yf, dy, xcte)
            implicit none
            real :: yo, yf, dy, suma
            ! función del integrando
            real :: fy, xcte
            integer :: N, i

            N = (yf-yo)/dy + 1
            suma = 0.0
            do i = 1, N-1
                suma = suma + fy(xcte, yo+i*dy)
            end do
            trapeciofy_y = (0.5*dy)*(fy(xcte,yo) + fy(xcte,yf) + 2*suma)
            return
        end function

        ! regla del trapecio para integrar funciones fz(cte,y)
        real function trapeciofz_y(yo, yf, dy, xcte)
            implicit none
            real :: yo, yf, dy, suma
            ! función del integrando
            real :: fz, xcte
            integer :: N, i

            N = (yf-yo)/dy + 1
            suma = 0.0
            do i = 1, N-1
                suma = suma + fz(xcte, yo+i*dy)
            end do
            trapeciofz_y = (0.5*dy)*(fz(xcte,yo) + fz(xcte,yf) + 2*suma)
            return
        end function