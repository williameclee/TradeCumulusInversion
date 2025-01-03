      subroutine sg1inv( pr, tr, nz, igamma, g, s, im, io, w, u, ierr )
      integer            nz, igamma, im, io, ierr
      SINGLE_OR_DOUBLE   pr, tr, g(0:nz), s(-1:nz+1), w(-1:nz+1,4)
      SINGLE_OR_DOUBLE   u(-1:nz+1)
c
c   Purpose
c
c       This routine solves the one-dimensional (X,Y-independent)
c       semigeostrophic invertibility relation for the Bernoulli 
c       function  Mstar  in geostrophic/isentropic coordinates.
c       All variables are nondimensional (for details, see  sg2inv).
c
c   On input
c
c       pr      Pressure ratio:  a = (p(bot) - p(top))/p(bot)
c
c       tr      Theta ratio:  b = Theta(bot)/(Theta(top) - Theta(bot))
c
c       nz      Number of grid intervals in  z
c
c       igamma  Specifies the approximation to  Gamma.  If  igamma=3,
c               Gamma  will be computed as part of the solution;
c               otherwise, the specified values in  g  will be used.
c
c       g       Array of size (0:nz) for values of  Gamma = d(PI)/dp.  
c               On input,  g  must contain the specified values of
c               Gamma  (unless  igamma=3,  in which case an initial
c               approximation will be computed internally).
c
c       s       Array of size (-1:nz+1) containing the input data:
c               s(nz+1) = PI(top)
c               s(   k) = sigmastar  at  z = k/nz,  k=0,...,nz
c               s(  -1) = phi(bot)
c
c       im      Specifies method (see description below):
c               1: trapezoidal       2: centered difference
c
c       io      Specifies output:  
c               io=0  no output
c               io=1  print convergence trace for Newton iteration
c
c       w       Work array of size at least (-1:nz+1,4)
c
c   On return
c
c       u       Contains values of  M:  u(k) = M  at  k/nz
c
c       g       Contains corresponding values of  Gamma = d(PI)/dp
c               (unchanged from input unless  igamma=3)
c
c       ierr    Error flag (zero on normal return)
c
c   Error condition
c
c       ierr > 0:  The solution did not converge in  ierr  iterations;
c                  u  and  g  may not be accurate (igamma=3 only).
c
c   Methods
c
c       Two methods are available for solving the problem as follows.
c
c       Trapezoidal (im=1):  First,  PI  is computed from  sigmastar
c       by solving downward in  z  with a trapezoidal discretization.
c       Then  Mstar  is computed from  PI  by solving upward in  z,
c       again with a trapezoidal discretization.
c
c       Centered difference (im=2):  The problem is discretized by 
c       second-order centered finite differences (using ghost points 
c       in the boundary conditions) and solved by Newton's method.  
c       An first approximation is computed by the trapezoidal method.
c
      integer          k, loop, nit
      SINGLE_OR_DOUBLE a, b, cnorm, dpi, fact1, fact2, fk, hz, hz2
      SINGLE_OR_DOUBLE kappa, mult, p, pi, tol, twohz, twohb, unorm
      parameter        ( tol = 1.0e-05, nit = 5 )
      logical          newgam
c
c   preliminaries
c
      ierr   = 0
      a      = pr
      b      = tr
      hz     = 1.0/nz
      hz2    = hz*hz
      twohz  = 2.0*hz
      twohb  = twohz/b
      kappa  = 2.0d0/7.0d0
      newgam = igamma.eq.3
c
c ************************* Trapezoidal Method *************************
c
c   compute  PI  in  w(.,1)  from  sigmastar
c
      fact1 =        kappa *a*hz/2.0
      fact2 = (1.0 - kappa)*a*hz/2.0
      pi = s(nz+1)
      w(nz,1) = pi
      if ( newgam ) g(nz) = 1.0/(pi*pi*sqrt( pi ) )
      do 30 k=nz-1,0,-1
          if ( newgam ) then
              g(k)= g(k+1)
                  if ( io.gt.0 ) write (7,2100) 0, g(k), pi
              do 10 loop=1,nit
                  fk = pi - w(k+1,1) - fact1*(g(k)*s(k)+g(k+1)*s(k+1))
                  dpi = fk/(1.0 + fact2*g(k)*s(k)/pi)
                  pi = pi - dpi
                  if ( pi.le.0.0 ) then
                      write (7,8000) pi
                      pi = 0.01
                  end if
                  g(k) = 1.0/(pi*pi*sqrt( pi ) )
                  if ( io.gt.0 ) write (7,2100) loop, g(k), pi, dpi
                  if ( abs( dpi ).le.tol*pi ) go to 20
   10         continue
              write (7,*) 'Newton iteration for initial  PI  failed'
   20         w(k,1) = pi
          else
              w(k,1) = w(k+1,1) + fact1*(g(k)*s(k) + g(k+1)*s(k+1))
          end if
   30 continue
c
c   compute  Mstar  in  u  from  PI
c
      fact1 = hz/(2.0*kappa*a)
      u(0) = s(-1) + b*w(0,1)/(kappa*a)
      do 40 k=1,nz
   40 u(k) = u(k-1) + fact1*(w(k-1,1) + w(k,1))
      u(  -1) = u(   1) + twohb*(s(-1) - u(0))
      u(nz+1) = u(nz-1) + twohz*s(nz+1) 
      if ( im.eq.1 ) return
c
c ********************* Centered Difference Method *********************
c
c   main loop for Newton iteration
c
      fact1 = 0.0
      do 100 loop=1,nit
c
c   set up the linear system for the correction
c
      w(-1,1) =  1.0/twohb
      w(-1,2) =  1.0
      w(-1,3) = -1.0/twohb
      w(-1,4) = s(-1) - (u(0) - (u(1) - u(-1))/twohb)
      do 50 k=0,nz
          if ( newgam ) then
              pi  = kappa*a*(u(k+1) - u(k-1))/twohz
              if ( pi.le.0.0 ) then
                  write (7,8000) pi
                  pi = 0.01
              end if
              g(k)  = 1.0/(pi*pi*sqrt( pi ) )
              p     = pi/g(k)
              fact1 = (kappa - 1.0)*a/(twohz*p)
          end if
          w(k,1) = -1.0/hz2 + fact1*s(k)
          w(k,2) =  2.0/hz2
          w(k,3) = -1.0/hz2 - fact1*s(k)
          w(k,4) = g(k)*s(k) + (u(k-1) - 2.0*u(k) + u(k+1))/hz2
   50 continue
      w(nz+1,1) = -1.0/twohz
      w(nz+1,2) =  0.0
      w(nz+1,3) =  1.0/twohz
      w(nz+1,4) = s(nz+1) - (u(nz+1) - u(nz-1))/twohz
c
c   eliminate the vertical ghost-point values from the system
c
      mult = w(0,1)/w(-1,1)
      w(0,2) = w(0,2) - mult*w(-1,2)
      w(0,3) = w(0,3) - mult*w(-1,3)
      w(0,4) = w(0,4) - mult*w(-1,4)
      mult = w(nz,3)/w(nz+1,3)
      w(nz,1) = w(nz,1) - mult*w(nz+1,1)
      w(nz,4) = w(nz,4) - mult*w(nz+1,4)
c
c   solve the system for the Newton correction by Gaussian elimination
c
      do 60 k=1,nz
          mult = w(k,1)/w(k-1,2)
          w(k,2) = w(k,2) - mult*w(k-1,3)
          w(k,4) = w(k,4) - mult*w(k-1,4)
   60 continue
      w(nz,4) = w(nz,4)/w(nz,2)
      do 70 k=nz-1,0,-1
   70 w(k,4) = (w(k,4) - w(k,3)*w(k+1,4))/w(k,2)
c
c   set the vertical ghost-point values of the correction
c
      w(-1,4) = (w(-1,4) - w(-1,2)*w(0,4) - w(-1,3)*w(1,4))/w(-1,1)
      w(nz+1,4) = (w(nz+1,4) - w(nz+1,1)*w(nz-1,4))/w(nz+1,3)
c
c   update the solution
c
      do 80 k=-1,nz+1
   80 u(k) = u(k) + w(k,4)
c
c   check for convergence
c
c     if ( .not.newgam ) return
      cnorm = 0.0
      unorm = 0.0
      do 90 k=-1,nz+1
      cnorm = max( cnorm, abs( w(k,4) ) )
      unorm = max( unorm, abs( u(k)   ) )
   90 continue
      if ( io.gt.0 ) write (7,2000) loop, cnorm, unorm
      if ( cnorm.le.tol*unorm ) return
  100 continue
      ierr = nit
      return
c
 2000 format('Newton iteration',i2,
     +       ':  cnorm =',1pe10.3,'  unorm =',e10.3)
 2100 format('Newton iteration',i2,
     +       ':  Gamma =',1pe10.3,'     PI =',e10.3:'  dPI =',e10.3)
 8000 format('sg1inv warning:  PI =',1pe10.2,'<=0:  reset to  0.01')
      end
