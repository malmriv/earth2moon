c     AUTHOR: Manuel Almagro (malmriv@correo.ugr.es). This is a preliminary
c     version to test my implementation of the RK4 method.

      program flymetothemoon
      implicit none
      real*8 phi, r_t, phi_t, pr_t, pr, pphi_t, p_phi, t
      real*8 dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      common /constants/ dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta

c     Define constants (these values are common to every function/subroutine)
      G = 6.67e-11
      MT = 5.9736e24
      ML = 0.07349e24
      dTL = 3.844e8
      omega = 2.6617e-6

c     Define and/or normalise necessary constants (these are common too)
      mu = ML/MT
      r = r/dTL
      r_tilde = dsqrt(1.d0+r**2.d0 - 2.d0*r*cos(phi-omega*t))
      delta = G*MT/dTL**3.d0

c     Define initial conditions. We can suppose the rocket at the surface of
c     the Earth,


      end program flymetothemoon

c     RK4 method
      function rk4(x0,y0,x,h) result(ans)
      implicit none
      real*8 pphi_t, y0, y, x0, x, h
      real*8 k1, k2, k3, k4, ans
      integer N, i
      N = nint((x-x0)/h)
      y = y0
      do i=1,N
        k1 = h*pphi_t(x0, y)
        k2 = h*pphi_t(x0+0.5*h, y+0.5*k1)
        k3 = h*pphi_t(x0+0.5*h, y+0.5*k2)
        k4 = h*pphi_t(x0+h, y+k3)
        y = y + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
        x0 = x0 + h
      end do
      ans = y
      end function rk4

c     Differential eqs. to solve.
      function pphi_t(phi,t) result(ans)
      implicit none
      real*8 phi, t, ans
      real*8 dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      common /constants/ dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      ans = (-delta*mu*r/r_tilde**3.d0)*sin(phi-omega*t)
      end function pphi_t

      function phi_t(pphi,t) result(ans)
      implicit none
      real*8 pphi, t, ans
      real*8 dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      common /constants/ dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      ans = pphi/r**2.d0
      end function phi_t

      function pr_t(pphi,phi,t) result(ans)
      implicit none
      real*8 pphi,phi, t, ans
      real*8 dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      common /constants/ dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      ans = pphi**2.d0/r**3.d0 -delta*(1/r**2.d0 + mu/r_tilde**3.d0*(
     &      r - cos(phi-omega*t)))
      end function pr_t


      function r_t(pr,t) result(ans)
      implicit none
      real*8 pr, t, ans
      real*8 dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      common /constants/ dTL, MT, ML, m, G, omega, mu, r, r_tilde, delta
      ans = pr
      end function r_t
