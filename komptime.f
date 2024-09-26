c This routine calculates time evolution of Kompaneets
c solution for a wind-blown bubble
c units are [L] = H, [rho] = rho_0, [L] = L_0
      program timedep
c
      real gamma
      real y,t,y0,t0,dy,dt,omega,omeganew,eth,ethnew,pressure
      open(1,file='result',status='unknown')

c
c  set constants
      gamma = 5./3.
c  set initial values for bubble from spherical similarity solution
      t0 = 3.38e-2
      y0 = 0.102
      omega = 4.189e-3
      eth  = 1.536e-2
      t = t0
      y = y0
c  set step size
      dy = 0.001
c  set stop condition
      yfinal = 1.98
c      yfinal = 1.999
c
  100 continue
      dt = dy/sqrt(0.5*(gamma**2-1.)*eth/omega)
c  for radiative shell, use following expression
c      dt = dy/sqrt((gamma-1.)*eth/omega)
      y = y + dy
      t = t + dt
c  get new volume
      omeganew = volume(y)
c  get new energy; first get old pressure, then 1st order time stepping
      pressure = (gamma - 1.)*eth/omega
      ethnew = eth + dt - pressure*(omeganew - omega)
c
      write(1,1000) y,t,ethnew,omeganew
 1000 format(5(1pe11.4))
c  reset old values to be the current new ones and start loop over
c  unless stop condition is met
      if (y.GE.yfinal) go to 200
      omega = omeganew
      eth=ethnew
      go to 100
c
  200 continue
c
      stop
      end
c
c
      real function radius(z,y)
c
      real z,y
      radius = 2.*acos(0.5*exp(0.5*z)*(1. - 0.25*y*y + exp(-z)))
      return
      end
c
      real function volume(y)
c
      integer npts
      real z1,z2,pi
c
      pi = 2.*asin(1.0)
c  get endpoints of integration z1 and z2
      z1 = -2.*alog(1. - 0.5*y)
      z2 = -2.*alog(1. + 0.5*y)
c  integrate over all cross sections
      npts=100
      dz = (z1-z2)/float(npts)
      volume = 0.
      do 10 n=1,npts
      zn = z2 + (2.*n-1.)*0.5*dz
      volume = volume + (radius(zn,y))**2
   10 continue
      volume = volume*pi*dz
c
      return
      end

