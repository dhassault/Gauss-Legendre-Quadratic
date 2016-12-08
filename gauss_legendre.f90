program gauss_legendre
    use functions
    implicit none

    character (len=14) :: cf
    integer :: k,n,npoly
    double precision:: a,b,gl,x
    double precision, allocatable :: xk(:), wk(:),f(:)

    write(6,*) 'This program is using the method of Guass-legendre.'

!**************************************************************

    !the bounded interval [a,b]:
    a=-1
    b=1

    !the file where the coefficients are stored:
    cf="cf.txt"
    npoly=2 !don't forget to change the cf.txt file also

    do n=1, 15
      allocate (xk(n), wk(n), f(n))
      call GL_xw (n, xk, wk)

      do k = 1, n
        f(k)=2/(1+poly(npoly,cf,xk(k)))
      end do

      call Quad_GL(f,a,b,n,xk,wk)
      deallocate (xk, wk, f)
    end do

!************************************************************

    !Generalization of the algorithm on [a,b]
    !for validation of the code

    !the new bounded interval is:
    a=-acos(-1.d0)
    b=acos(-1.d0)

    do n=1, 15
      allocate (xk(n), wk(n), f(n))
      call GL_xw (n, xk, wk)

      do k = 1, n
        f(k)=sin(((b-a)*0.5)*xk(k)+((a+b)*0.5))*sin(((b-a)*0.5)*xk(k)+((a+b)*0.5))
      end do


      gl=0
      do k = 1, n
        if ( xk(k)>a .AND. xk(k)<b)  gl=gl+wk(k)*f(k)
      end do
      gl=((b-a)*0.5)*gl
      write(6,*) '        Ordre     Result'
      write(6,*) k,    gl

      deallocate (xk, wk, f)
    end do

!***************************************************************

end program gauss_legendre
