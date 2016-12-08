module functions
  implicit none

contains
!*********************************************************

  subroutine GL_xw(n, xk, wk)
    !Calculating the roots xk(n) with Newton-Raphson
    integer :: i,n,k
    double precision:: pi,eps,P0,P1,P,dP,x,dx
    double precision:: xk(n), wk(n)
    pi=acos(-1.d0)
    eps=1e-10
    do k=1,n
        x=cos((4.d0*k-1.d0)*pi/(4.d0*n+2.d0))
        do
            P0 = 1.d0
            P1 = x
            do i = 2, n
                P=(((2.0d0*i-1.0d0)*x*P1-(i-1.0d0)*P0))/i
                P0 = P1
                P1 = P
            end do
            dP=(n*P0-n*x*P1)/(1.d0-x*x)
            dx = P/dP
            x=x-dx
            if (abs(dx)<eps) exit
        end do
        xk(k) = x
        wk(k)=2.d0/((1.d0-x*x)*(dP*dP))
    end do
  end subroutine GL_xw

!*****************************************************

  subroutine Quad_GL(f,a,b,n,xk,wk)
    !Evaluating the quadratic Gauss-Legendre itself
    integer :: n,k
    double precision:: a,b,gl
    double precision:: f(n), xk(n), wk(n)

    gl=0
    do k = 1, n
      if ( xk(k)>a .AND. xk(k)<b)  gl=gl+wk(k)*f(k)
    end do
    !write(6,*) gl

  end subroutine Quad_GL

!*******************************************************

  function poly(npoly,cf,x)
    !Evaluating the polynomial at npoly order with Horner Method
    !coefficients of the polynomial are in cf.txt
    character (len=14) :: cf
    logical :: existing
    integer :: i,npoly
    double precision :: a,b,poly,x

    inquire(file=cf,exist=existing)
    if    ( .NOT. existing ) then
      write(*,*)"Cannot open the file..."
      STOP
    endif

    open(10, file=cf, action="read")
    read(10,*) a,b

    poly=a*x+b
    do i=npoly-2,0,-1
      read(10,*) a
      poly=poly*x+a
    end do
    close(10)
  end function poly

!***********************************************************
end module functions
