c     AUTHOR: Manuel Almagro (malmriv@correo.ugr.es). This is a preliminary
c     version to test my implementation of the RK4 method.

      program flymetothemoon
      implicit none
      real*8 y(1:4), yaux(1:4), k1(1:4), k2(1:4), k3(1:4), k4(1:4)
      real*8 p_phi, phi, p_r
      real*8 dTL,MT,RT,ML,m,G,omega,mu,delta
      real*8 theta, t, t_total, h
      integer i
      common /constants/ dTL,MT,RT,ML,m,G,omega,mu,delta

c     Define parameters (these values are common to every function/subroutine)
      G = 6.67d-11
      MT = 5.9736d24
      ML = 7.349d22
      RT = 6.378160d6
      dTL = 3.844d8
      omega = 2.6617d-6
      m = 4.5702d4 !Apollo 11 launch mass
      theta = 0.d0 !Polar angle of the Moon
      h = 1d-2
      t = 0.d0
      t_total = 1800

c     Define initial conditions.
      y(1) = RT !Initial r
      y(2) = 0.d0 !Initial phi, 30deg
      y(3) = 0.9*dsqrt(2.d0*G*MT/RT)*cos(theta-y(2))/dTL !Initial p_r
      y(4) = (RT/dTL)*dsqrt(2.d0*G*MT/RT)/dTL*sin(theta-y(2)) !Initial p_phi

c     Define and/or normalise necessary constants (these are common, too)
      mu = ML/MT
      delta = G*MT/(dTL**3.d0)
      y(1) = y(1)/dTL

      open(20,file="./trajectory.txt")

      do i=1,nint(t_total/h)

c     Solve for r with RK4:
      k1(1) = h*y(3)
      k2(1) = h*(y(3)+k1(1)/2.d0)
      k3(1) = h*(y(3)+k2(1)/2.d0)
      k4(1) = h*(y(3)+k3(1))
      yaux(1) = y(1) + 1.d0/6.d0*(k1(1)+2.d0*k2(1)+2.d0*k3(1)+k4(1))

c     Solve for phi with RK4:
      k1(2) = h*phi(y(4),y(1),t)
      k2(2) = h*phi(y(4)+k1(2)/2.d0,y(1),t+h/2.d0)
      k3(2) = h*phi(y(4)+k2(2)/2.d0,y(1),t+h/2.d0)
      k4(2) = h*phi(y(4)+k3(2),y(1),t+h)
      yaux(2) = y(2) + 1.d0/6.d0*(k1(2)+2.d0*k2(2)+2.d0*k3(2)+k4(2))

c     Solve for p_r with RK4:
      k1(3) = h*p_r(y(1),t,y(4),y(2))
      k2(3) = h*p_r(y(1)+k1(3)/2.d0,t+h/2.d0,y(4),y(2))
      k3(3) = h*p_r(y(1)+k2(3)/2.d0,t+h/2.d0,y(4),y(2))
      k4(3) = h*p_r(y(1)+k3(3),t+h,y(4),y(2))
      yaux(3) = y(3) + 1.d0/6.d0*(k1(3)+2.d0*k2(3)+2.d0*k3(3)+k4(3))

c     Solve for p_phi with RK4:
      k1(4) = h*p_phi(y(2),y(1),t)
      k2(4) = h*p_phi(y(2)+k1(4)/2.d0,y(1),t+h/2.d0)
      k3(4) = h*p_phi(y(2)+k2(4)/2.d0,y(1),t+h/2.d0)
      k4(4) = h*p_phi(y(2)+k3(4),y(1),t+h)
      yaux(4) = y(4) + 1.d0/6.d0*(k1(4)+2.d0*k2(4)+2.d0*k3(4)+k4(4))

c     Save results, update values, update time
      y(1) = yaux(1)
      y(2) = mod(yaux(2),2.d0*3.14159265359)
      y(3) = yaux(3)
      y(4) = yaux(4)
      write(20,*) y(1), y(2), omega*t
      t = t+h
      end do
      close(20)
      end program flymetothemoon


c     Differential eqs. to solve.
      function p_phi(phi,r,t) result(ans)
      implicit none
      real*8 phi, r, t, ans, r_tilde
      real*8 dTL,MT,RT,ML,m,G,omega,mu,delta
      common /constants/ dTL,MT,RT,ML,m,G,omega,mu,delta
      r_tilde = dsqrt(1.d0+r**2.d0 - 2.d0*r*cos(phi-omega*t))
      ans = -delta*mu*r/r_tilde**3.d0*sin(phi-omega*t)
      end function p_phi

      function phi(pphi,r,t) result(ans)
      implicit none
      real*8 pphi, r, t, ans
      real*8 dTL,MT,RT,ML,m,G,omega,mu,delta
      common /constants/ dTL,MT,RT,ML,m,G,omega,mu,delta
      ans = pphi/(r**2.d0)
      end function phi

      function p_r(r,t,pphi,phi) result(ans)
      implicit none
      real*8 r, t, pphi, phi, ans, r_tilde
      real*8 dTL,MT,RT,ML,m,G,omega,mu,delta
      common /constants/ dTL,MT,RT,ML,m,G,omega,mu,delta
      r_tilde = dsqrt(1.d0+r**2.d0 - 2.d0*r*cos(phi-omega*t))
      ans = pphi**2.d0/r**3.d0 -delta*(1/r**2.d0 + mu/r_tilde**3.d0*(
     &      r - cos(phi-omega*t)))
      end function p_r
