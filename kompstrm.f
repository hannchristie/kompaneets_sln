c  This program calculates the streamlines for the
c  Kompaneets solution 
c
      program strmline 
c
      integer ny,npts
      parameter(npts=100,ny=189)
      real y0,y,dy,yold 
      real r(npts,ny),z(npts,ny),slope(npts,ny),theta(npts)
c  naglib related declarations
      integer ifail,nfmax
      real*8 x,eps,eta,f,a,b
      real*8 rnag,znag,snag,ynag
      common/com1/rnag,znag,snag,ynag
      external f
      open(1,file='rstrm',status='unknown')
      open(2,file='zstrm',status='unknown')
c
c  set y0, dy 
      y0 = 0.10
      dy = 0.01
      y=y0
c
c  at initial y0, find z and r of starting points
c  pick initial z's by dividing initial curve into vertical strips
      z1 = -2.*alog(1. - 0.5*y)
      z2 = -2.*alog(1. + 0.5*y)
      dz = (z1-z2)/float(npts-9)
      z(1,1) = z2 + dz
      do 10 n=2,npts-10 
      z(n,1) = z(n-1,1) + dz
   10 continue
c  last 10 points very close to axis
      do 15 n=npts-9,npts-1
      z(n,1) = z1*(1. - (npts+1-n)*0.0005)
      write(6,*) z1,z(npts-10,1),z(n,1)
   15 continue
      z(npts,1) = z1*(1. - 0.0001)
      write(6,*) z1,z(npts-10,1),z(npts,1)
      do 20 n=1,npts
      r(n,1) = radius(z(n,1),y)
   20 continue
      do 25 n=npts-10,npts
      theta(n) = atan(r(n,1)/z(n,1))*180./3.14159
      write(6,*) n,'theta',theta(n)
   25 continue
c  originate straight line at each point and find new z for ynew = y + dy 
      do 30 m = 2,ny
      yold = y
      y = y + dy
      do 40 n = 1,npts
      slope(n,m-1) = -1./drdz(z(n,m-1),yold)
c  use nonlinear solver to find zsol such that
c  r(n,m-1) + slope(n,m-1)*zsol = radius(zsol,y), where y = y+dy
c  and radius = 2.*acos(0.5*exp(z/2.)*(1 - 0.25*y*y + exp(-z))
      x = z(n,m-1)
      eps = 1.e-9
      eta = 0. 
      ifail = 0
      nfmax = 1000
c
      rnag = dble(r(n,m-1))
      znag = x
      snag = dble(slope(n,m-1))
      ynag = dble(y)
c
      call c05ajf(x,eps,eta,f,nfmax,ifail)
c
      z(n,m) = x
      r(n,m) = radius(z(n,m),y)
c
   40 continue      
   30 continue
c  plot r(n,*) vs. z(n,*) to get streamline in r,z plane for nth point on 
c  initial sphere
c  print out values
      do 60 m=1,ny
      do 70 n=1,npts
      write(1,1000) r(n,m)
   70 continue
   60 continue
      do 80 m=1,ny
      do 90 n=1,npts
      write(2,1000) z(n,m)
   90 continue
   80 continue
 1000 format(1pe11.4)

      
      stop
      end


      real function radius (z,y)
c
      real z,y
c
      radius = 2.*acos(0.5*exp(z/2.)*(1. - 0.25*y*y + exp(-z)))
c
      return
      end
c
      real function drdz (z,y)
c
      real z,y,rad
c
      rad = 2.*acos(0.5*exp(z/2.)*(1. - 0.25*y*y + exp(-z)))
      drdz = (exp(-z/2) - cos(rad/2.))/sin(rad/2.)
c
      return
      end
c
      real*8 function f (x)
c
      real*8 x
      intrinsic exp
      intrinsic cos
      real*8 rnag,znag,snag,ynag
      common/com1/rnag,znag,snag,ynag

c
      f = cos(0.5*(rnag + snag*(x-znag))) 
     #- 0.5*exp(x/2)*(1. - 0.25*ynag*ynag + exp(-x))
c
      return
      end
