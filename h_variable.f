      program flymetothemoon
      implicit none
      real*8 y(1:4),f(1:4),k1(1:4),k2(1:4),k3(1:4),k4(1:4)
      real*8 y_test(1:4), y_aux(1:4)
      real*8 deg2rad, v_e, h, h_max, err, err_max, s
      real*8 dTL,MT,RT,ML,G,omega,mu,delta
      real*8 theta, t, t_total, hamilt, r_prima
      integer i
      common /constants/ dTL,MT,RT,ML,G,omega,mu,delta

c     Define parameters (these values are common to every function/subroutine)
      G = 6.67e-11
      MT = 5.972e24
      ML = 7.3477e22
      RT = 6.378160e6
      dTL = 3.844e8
      omega = 2.6617e-6 !Angular frequency of the Moon

c     These parameters control the simulation
      theta = deg2rad(59.0d0)
      v_e = dsqrt(2.d0*G*MT*(1/RT-1/dTL))/dTL
      h = 1.d0
      err_max = h**5.d0
      t = 0.d0
      t_total = 7.d0*86400.d0 !A week

c     Define initial conditions.
      y(1) = RT/dTL !r
      y(2) = deg2rad(53.d0) !phi
      y(3) = v_e*cos(theta-y(2)) !dr/dt
      y(4) = y(1)*v_e*sin(theta-y(2)) !d(phi)/dt
      f = 0.d0

c     Define and/or normalise necessary constants (these are common, too)
      mu = ML/MT
      delta = G*MT/(dTL**3.d0)

      open(20,file="./trajectory.txt")
      open(30,file="./hamilt.txt")
      do while(t .lt. t_total)
c       Runge-Kutta (4th order) method. We test results with h and with h/2
        call functions(y,t,f)
        k1 = h*f
        call functions(y+k1/2.d0,t+h/2.d0,f)
        k2 = h*f
        call functions(y+k2/2.d0,t+h/2.d0,f)
        k3 = h*f
        call functions(y+k3,t+h,f)
        k4 = h*f
        y_aux = y + 1.d0/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)

        call functions(y,t,f)
        k1 = (h/2.d0)*f
        call functions(y+k1/2.d0,t+h/4.d0,f)
        k2 = (h/2.d0)*f
        call functions(y+k2/2.d0,t+h/4.d0,f)
        k3 = (h/2.d0)*f
        call functions(y+k3,t+(h/2.d0),f)
        k4 = (h/2.d0)*f
        y_test = y + 1.d0/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)

c       This is where the decision to change the value of h is taken
55      err = (16.d0/15.d0)*maxval(abs(y_test-y_aux))
        s = max((err/err_max)**0.2d0,1d-8)
        h_max = h/s
        if(s .gt. 2.d0) then
          h = h_max
          goto 55
        end if
        if(s .lt. 2.d0) then
          h = h/2.d0
          t = t + h
          y = y_test
        end if
        if(h .lt. h_max) h = 2.d0*h

c       Save results
        write(20,*) y, omega*t, t

c       Check for H conservation
        r_prima = dsqrt(1.d0+y(1)**2.d0-2.d0*y(1)*cos(y(2)-omega*t))
        hamilt = ((y(3)*dTL)**2.d0)/2.d0 + ((y(4)*dTL)**2.d0)/
     &  (2.d0*y(1)**2.d0)-G*MT/(y(1)*dTL)-G*ML/(r_prima*dTL)
        write(30,*) hamilt
      end do
      close(20)
      close(30)
      end program flymetothemoon


c     Differential eqs. to solve.

      subroutine functions(y, t, f)
      implicit none
      real*8 f(1:4), y(1:4), r_prima, t
      real*8 dTL,MT,RT,ML,G,omega,mu,delta
      common /constants/ dTL,MT,RT,ML,G,omega,mu,delta

c     Rocket-Moon distance (changes constantly unlike mu and delta)
      r_prima = dsqrt(1.d0+y(1)**2.d0-2.d0*y(1)*cos(y(2)-omega*t))

      f(1) = y(3)
      f(2) = y(4)/y(1)**2.d0
      f(3) = (y(4)**2.d0)/(y(1)**3.d0) - delta*(1.d0/(y(1)**2.d0)+mu/
     &(r_prima**3.d0)*(y(1)-cos(y(2)-omega*t)))
      f(4) = -(delta*mu*y(1)/(r_prima**3.d0))*sin(y(2)-omega*t)

      return
      end subroutine functions

c     Some functions to make testing easier
      function deg2rad(deg) result(rad)
      implicit none
      real*8 deg, rad
      rad = deg*(3.1415926535897932/180.d0)
      end function deg2rad
