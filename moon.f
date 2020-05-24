c     AUTHOR: Manuel Almagro (malmriv@correo.ugr.es). This is a preliminary
c     version to test my implementation of the RK4 method.

      program flymetothemoon
      implicit none
      real*8 x0, y0, x, h, rk4
c     Establish initial conditions
      x0 = 1.d0
      y0 = 3.d0
c     Point to evaluate and step
      x = 2.d0
      h = 1E-4
c     Call RK4 with these parameters
      write(*,*) rk4(x0,y0,x,h)
      end program flymetothemoon

c     RK4 method
      function rk4(x0,y0,x,h) result(ans)
      implicit none
      real*8 func, y0, y, x0, x, h
      real*8 k1, k2, k3, k4, ans
      integer N, i
      N = nint((x-x0)/h)
      y = y0
      do i=1,N
        k1 = h*func(x0, y)
        k2 = h*func(x0+0.5*h, y+0.5*k1)
        k3 = h*func(x0+0.5*h, y+0.5*k2)
        k4 = h*func(x0+h, y+k3)
        y = y + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
        x0 = x0 + h
      end do
      ans = y
      end function rk4

c     The differential eq. to solve is dy/dx = x+y/x
c     With y(1) = 3, analitical solution is y = x(x+2), y(2) = 8
      function func(x,y) result(ans)
        implicit none
        real*8 x, y, ans
        ans = x+y/x
      end function func
