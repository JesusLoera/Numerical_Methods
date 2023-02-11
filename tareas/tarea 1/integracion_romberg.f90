! Programa elaborado por Jesús Eduardo Loera Casas
! Fecha 01/02/23

! Aquí ingresa tu función.
        REAL FUNCTION fx(x)
            IMPLICIT NONE
            REAL :: x
            ! Parametros
            REAL :: mo, g, q, u
            g = 9.8; u = 1800; mo = 160000; q = 2500
            ! ecuación del cohete
            fx = u*log((mo)/(mo-q*x))-g*x
            RETURN 
        END

        REAL FUNCTION romberg_00(a, b)   
            IMPLICIT NONE
            REAL :: a, b, fx
            romberg_00 = 0.5*(b-a)*(fx(a)+fx(b))
            RETURN
        END

        RECURSIVE FUNCTION romberg_n0(n, a, b) result(integral)
            IMPLICIT NONE
            INTEGER :: n, k
            REAL :: a, b, fx, romberg_00, hn, sum, integral
            IF (n==0) THEN
                integral = romberg_00(a, b)
            ELSE 
                hn = (b-a)/(2**n)
                sum = 0
                DO k = 1, 2**(n-1)
                    sum = sum + fx(a + (2*k-1)*hn)
                END DO
                integral = 0.5*(romberg_n0(n-1, a, b)) + hn*sum
            END IF
            RETURN
        END

        REAL FUNCTION romberg_nm(n, m, a, b)
            IMPLICIT NONE
            INTEGER :: n, m, i, j
            REAL :: a, b, romberg_ij, aux, romberg_00, romberg_n0, tol
            REAL, DIMENSION ((n+1),(m+1)) :: matrix

            ! Definimos una tolarancia al error con cada iteración
            tol = 0.00001
            
            matrix(1,1) = romberg_00(a,b)

            DO i = 2, n + 1
                matrix(i,1) = romberg_n0(i-1, a, b)
                ! Evaluamos el criterio de convergencia
                IF (abs(matrix(i,1)-matrix(i-1,1)) .le. tol ) THEN
                    write(*,*) "La integral convergio con una &
     &tolerancia de: ", tol
                    romberg_nm = matrix(i,1)
                    RETURN
                END IF
            END DO

            DO j = 2, m+1
                DO i = j, n+1
                    aux = (4**(j-1))*matrix(i, j-1) - matrix(i-1, j-1)
                    romberg_ij = (1.0/((4.0**(j-1))-1))*(aux)
                    matrix(i,j) = romberg_ij
                    ! Evaluamos el criterio de convergencia
                    IF (abs(matrix(i,j)-matrix(i-1,j)) .le. tol ) THEN
                        write(*,*) "La integral convergio."
                        romberg_nm = matrix(i,j)
                        RETURN
                    END IF
                END DO
            END DO

            write(*,*)"La integral no convergió con la tolerancia dada."
            romberg_nm = matrix(n+1, m+1)
            RETURN
        END

        PROGRAM main
 
            IMPLICIT NONE
            INTEGER :: n, m
            REAL :: a, b
            REAL :: romberg_nm, int_nm

            n = 10 ; m = 4
            a = 0  ; b = 30

            int_nm = romberg_nm(n, m, a, b)

            WRITE(*,"(A27,I2,A1,I2,A6,F20.10)") "La aproximacion de &
     &romberg(",n,",",m, ") es: ",  int_nm 

        END PROGRAM main