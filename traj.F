      program main 
c
c
c     This program uses a newton iteration to solve the invertibility
c     equation for the shallow water version of the hadley cell model.
c
      integer  jc, jt, kc, kcp, kcm, jl, is 
      parameter (kc = 40,  kcp=kc+1, kcm=kc-1)
      parameter (jc = 81, jt=2*jc-1, jl=60)
      parameter (is=19)
      integer  i, jm, jp, j, k, iter, nmax, n, ks, kp, km, itr
      integer  nn, inc, jump, js, mm, nup, niter, ierr, itm
      integer  shr, ehr, dhr
      SINGLE_OR_DOUBLE  x(-jc+1:jc-1,-1:kcp), l(-jc+1:jc-1,0:kc)
      SINGLE_OR_DOUBLE  d(-jc+1:jc-1,0:kc), u(-jc+1:jc-1,0:kc)
      SINGLE_OR_DOUBLE  g(-jc+1:jc-1,0:kc), phid(-jc:jc)
      SINGLE_OR_DOUBLE  sc(-jc:jc), scs(-jc:jc)
      SINGLE_OR_DOUBLE  fac1(-jc:jc), f1(-jc:jc)
      SINGLE_OR_DOUBLE  z(0:kc), theta(0:kc)
      SINGLE_OR_DOUBLE  gamma(-jc:jc,0:kc), sigma(-jc:jc,0:kc)
      SINGLE_OR_DOUBLE  q0n, qalfn, phicn, sinpcn, xerfn, degn
      SINGLE_OR_DOUBLE  q0s, qalfs, phics, sinpcs, xerfs, degs
      SINGLE_OR_DOUBLE  cs, eps, ae, thetat, thetab, pt, pb
      SINGLE_OR_DOUBLE  r, gr, omega, rom, cp, kappa
      SINGLE_OR_DOUBLE  pi, sqpi, dz, tdz, dt, alpha, beta
      SINGLE_OR_DOUBLE  sigma0, presb, gammab, pieb, msb, mst
      SINGLE_OR_DOUBLE  convd, convr, times
      SINGLE_OR_DOUBLE  eps1, zeta, zeta0, qhat, argmnt, tau
      SINGLE_OR_DOUBLE  xjs, fj, dj, dsc, rn, erf
      SINGLE_OR_DOUBLE  dxds, dxdz2, f2, f3, f4, term2
      SINGLE_OR_DOUBLE  xold, fac2, gkp, gkpp, xs, xsp, xsm, term3
      SINGLE_OR_DOUBLE  gkm, gkmp, dmz, pie, phimax, phimin
      SINGLE_OR_DOUBLE  t1, t2, ang, sigs, lambda
      SINGLE_OR_DOUBLE  zfac, lz
      real              csval, fp(is,4), tht(is), pr(is), temc(is)
      real              cot, thetak, pmin, pmax, tmin, tmax, lat
      real              qht(-jl:jl), latr, thet(0:kc)
      real              dl, thdot(-jc:jc,0:kc)
c     data (pr(i),i=1,is) /1000., 950., 900., 850., 800., 750.,
c    2     700., 650., 600., 550., 500., 450., 400., 350., 300., 250.,
c    3     200., 150., 125., 100./
c     data ( temc(i),i=1,is) /24.60, 21.28, 18.42, 15.96, 13.54,
c    2     10.78, 7.55, 3.98, 0.39, -3.28, -7.34, -12.08, -17.89,
c    3    -24.85, -33.33, -43.76, -56.00, -69.07, -73.76, -78.46/
      data (pr(i),i=1,is) /1000., 950., 900., 850., 800., 750.,
     2     700., 650., 600., 550., 500., 450., 400., 350., 300., 250.,
     3     200., 150., 100./
      data ( temc(i),i=1,is) /299.73, 301.73, 303.80, 306.02, 308.46,
     2     311.10, 313.90, 316.87, 320.21, 324.08, 328.09, 331.94,
     3     336.09, 340.81, 344.8, 347.2, 350.35, 356.91, 367.74/
      character gam*(1), sa*(1), insig*(1)
      common /param/ thet, qht(-jl:jl), lat(-jl:jl), thetat,
     2               thetab, pi
c
c   statement function for cot
c
      cot(ang) = 1. / tan(ang)
c
c   open file for output
c
      open(7, file='invert.out', status='unknown', form='formatted')
c
c   open workstation
c
      call opngks
c 
c   turn off clipping
c
      call gsclip( 0 )
C
C  SOME CONSTANTS
C
      R     = 287.d0
      CP    = 1004.d0
      KAPPA = R/CP
      gr    = 9.81d0
      AE    = 6.371d+06 
      OMEGA = 7.292d-05
      ROM   = AE * OMEGA
      pi    = acos(-1.d0)
      convd = 180.d0 / pi
      convr = pi / 180.d0
      sqpi  = sqrt(pi) 
      thetat= 360.d0
      thetab= 300.d0
      dt    = (thetat - thetab) / kc
      pb    = 100000.0d0
      pt    =  12500.0d0
      alpha = (pb-pt) / pb
      beta  = thetab / (thetat-thetab)
      sigma0= (pb-pt) / (thetat-thetab)
      do 1  i=1,is
         fp(i,1) = pr(i)*100.
         tht(i) = (temc(i)+273.15) * (pr(1) / pr(i))**kappa
         print * , i, tht(i), fp(i,1)
    1 continue
c
c   fit cubic spline to data
c
      call csset(is, tht, fp, 4, 4, ierr)
c
      cs    = alpha*r*(thetat-thetab)
      eps   = 4.d0 * rom*rom / cs
      print 2, eps, sqrt(cs)
    2 format('epsilon = ',1p,e10.3,' c = ',e10.3)
c
c   specify form for initial sigma
c
      print 3
    3 format(/'Enter form for initial sigma:',
     2       /'      c    for constant with height',
     3       /'      v    for variable with height ...>')
      read(5,*) insig
c
c   specify form of gamma
c
      print 5
    5 format(/'Enter form of gamma for problem:',
     2       /'      a    for anelastic case',
     3       /'      b    for boussinesq case',
     4       /'      g    for general case    ...>')
      read(5,*) gam
    6 print 7
    7 format(/'Enter form of sigma anomalies for problem:',
     2       /'      i    for internal anomalies ',
     3       /'      b    for boundary anomalies ...>')
      read(5,*) sa
      if(sa .ne. 'i' .and. sa .ne. 'b') then
         print 8, sa
         go to 6
      endif
    8 format(' Improper specification for sigma anamalies: ',a1)
C
C  HEATING PARAMETERS
C
c   northern heat source
c
      PRINT *,'Specify magnitude (K/day) of northern heating ...>'
      READ  *, q0n
      if(q0n .ne. 0.0) then
      PRINT *,'Specify central lat. (in deg) of northern heating ...>'
      READ  *, degn
       PRINT *,'Specify width (qalfn) of northern heating '
       PRINT *,'(e.g., qalfn=15, half-width of 8 degrees)'
       PRINT *,'(e.g., qalfn=30, half-width of 4 degrees)'
       READ  *, qalfn
      else
       degn = 0.0
       qalfn = 0.0
      endif
      phicn = degn / convd
      sinpcn = sin(phicn)
      xerfn = erf(qalfn*(1.+sinpcn)) + erf(qalfn*(1.-sinpcn))
c
c   southern heat source
c
      PRINT *,'Specify magnitude (K/day) of southern heating ...>'
      READ  *, q0s
      if(q0s .ne. 0.0) then
       PRINT *,'Specify central lat. (in deg) of southern heating ...>'
       READ  *, degs
       PRINT *,'Specify width (qalfs) of southern heating ...>'
       READ  *, qalfs
      else
       degs = 0.0
       qalfs = 0.0
      endif
      phics = degs / convd
      sinpcs = sin(phics)
      xerfs = erf(qalfs*(1.+sinpcs)) + erf(qalfs*(1.-sinpcs))
      print 9
    9 format('Specify lambda factor to control level of heating'
     2      ,/' lambda = 1.0, top of heating at 360K'
     3      ,/' lambda = 1.2, top of heating at 350K'
     4      ,/'                          your choice ...>')
      read *, lambda
C
C     ******************* CREATE COORDINATES ********************
C
C     sc(j) goes from -1 to 1
C
      dsc = 1.d0 / jc
      do 10 j=-jc,jc
         sc(j)   = j*dsc
         if(sc(j) .gt. 1.) sc(j) =  1.
         if(sc(j) .lt.-1.) sc(j) = -1.
         scs(j)  = sc(j)*sc(j)
         fac1(j) = eps*dsc*sc(j)
         f1(j)   = fac1(j)*(1.-scs(j))
         phid(j) = asin(sc(j))*convd
c        write(7,*) j, sc(j), phid(j)
   10 continue
      dz = 1.d0 / kc
      tdz = 2.d0*dz
      fac2= 2.d0*dsc*dz*dz
      do 20  k=0,kc
         z(k)     = k*dz
         theta(k) = thetab + k*dt
         thet(k)  = theta(k)      
   20 continue
c
c     ****************** INITIALIZE VARIABLES *******************
c
c   initialize pressure, gamma and x (i.e., m) at even points
c
      do 30  k=0,kc
         presb  = 1. - alpha*z(k)
         pieb   = presb**kappa
         if(gam .eq. 'a') then
            gammab = 1.d0
         else
            gammab = pieb**((kappa-1.)/kappa)
         endif
         msb    = (beta + (1.-presb**(1.+kappa))  
     2          / (alpha*(1.+kappa)) ) / (kappa*alpha)
         write(7,24)  k, presb, pieb, gammab, msb
   24    format(i4,1p,4e12.4)
         do 25 j=-jc+1,jc-1,2
            gamma(j,k) = gammab
            x(j,k)     = msb
   25    continue
   30 continue
      msb = x(0,1) - tdz*x(0,0)/beta
      mst = x(0,kcm) + tdz*(1.d0-alpha)**kappa / (kappa*alpha)
      do 32 j=-jc+1,jc-1,2
         x(j,-1) = msb
         x(j,kcp) = mst
   32 continue
c
c   initialize x (i.e., sin(phil)) at odd points
c
      do 40  k=-1,kcp
      do 40  j=-jc+2,jc-2,2
         x(j,k) = sc(j)
   40 continue
C
C     +*+*+*+*+*+*+*+*+*+ READ INPUT +*+*+*+*+*+*+*+*+*+*+*+*
c
c     specify timed (time in days that solution is to be computed for)
c     specify nmax (the number of relaxation sweeps to be performed)
c
      print *, 'Specify southern-most latitude of plots ...> '
      read *, phimin
      print *, 'Specify northern-most latitude of plots ...> '
      read *, phimax
   45 continue
      print 46
   46 format('Specify starting hour to compute trajectories for ',/,
     2        '(negative value to quit) ...>')
      read  *, shr 
      if(shr .lt. 0) go to 1000
      print *,'Specify ending hour to compute trajectories for ...>'
      read  *, ehr 
      print *,'Specify hour increment to compute trajectories at ...>'
      read  *, dhr
      print *,'Specify number of sweeps at each time ...>'
      read  *, nmax
      print *,'Specify how frequent to print residual ...>'
      read  *, nup
c     print *,'Specify number of newton iterations to do ...>'
c     read  *, niter 
      niter = 1
C
      print 51, shr, ehr, dhr
      print 52, degn, qalfn
      print 53, degs, qalfs
   51 format(' start hour  = ',i5', end hour = ',i5,
     2       ' hour increment = ',i5)
   52 format(' center of northern heating = ',f5.1,
     2       '  heating param. = ',f5.1)
   53 format(' center of southern heating = ',f5.1,
     2       '  heating param. = ',f5.1)
c
c   interpolate fields to actual latitude defined below
c
      dl = (phimax - phimin) / (2*jl)
      do 54  j=-jl,jl
         lat(j) = phimin + (j+jl)*dl
         latr = lat(j)*convr
         qht(j) = q0n*4.*qalfn*exp(-(qalfn*(sin(latr)-sinpcn))**2)
     2             / (xerfn*sqpi)
     3          + q0s*4.*qalfs*exp(-(qalfs*(sin(latr)-sinpcs))**2)
     4             / (xerfs*sqpi)
   54 continue
c
c   loop from shr to ehr
c
      do 999  itm=shr,ehr,dhr
      print 66, itm
   66 format(/,' computing the solution for hour ',i3,/)
      times = float(itm) * 3600.
c
c     ******************* COMPUTE SIGMA ***********************
c
      eps1 = 1.E-10
      DO 70 k=0,kc
         lz   = lambda * z(k)
         DO 60  j=-jc,jc
            qhat = q0n*4.*qalfn*exp(-(qalfn*(sc(j)-sinpcn))**2)
     2             / (xerfn*sqpi*86400.)
     3           + q0s*4.*qalfs*exp(-(qalfs*(sc(j)-sinpcs))**2)
     4             / (xerfs*sqpi*86400.)
            tau  = (qhat*times)/(thetat-thetab)
            argmnt = exp(-pi*tau)*tan(zeta/2.0)
         if(sa .eq. 'b') then
c  heating profiles producing anomalies at the boundaries
            if (lz  .lt. eps1) then
               zeta0 = 0.0
               sigma(j,k)= exp(-pi*tau)
               thdot(j,k) = 0.0	
            elseif (abs(zeta-pi) .lt. eps1) then
               zeta0 = pi
               sigma(j,k)= exp(pi*tau)
               thdot(j,k) = 0.0	
            else
               zeta0 = 2.*atan(argmnt)
               sigma(j,k)=sin(zeta0)/sin(zeta)
               thdot(j,k) = qhat*sin(pi*lz)*86400.
            end if
         else
c  heating profiles producing only internal anomalies
            if (lz .lt. eps1) then
               zeta0 = 0.0
               sigma(j,k)= 1.0
               thdot(j,k) = 0.0	
            elseif (lz .ge. 1.0) then
               zeta0 = z(k)
               sigma(j,k)= 1.0
               thdot(j,k) = 0.0	
            else
               zfac = 0.5*pi - atan(pi*lambda*tau+cot(pi*lz))
               sigma(j,k)= sin(zfac)**2 / sin(pi*lz)**2
               zeta0 = zfac / (pi * lambda)
               thdot(j,k) = qhat*sin(pi*lz)**2 * 86400.
            end if
         endif
         if(insig .eq. 'v') then
            thetak = (thetat - thetab)*zeta0 + thetab
            sigs = -csval(thetak, 1, is, thet, fp, 0) / sigma0
            sigma(j,k)= sigs * sigma(j,k)
         endif 
   60    continue
   70 continue
      if(itm .eq. 0) go to 290
C
C     ******************* BEGIN ITERATION *********************
C
      do 280  n=1,nmax
      rn = 0.0
c
c   zebra relaxation in horizontal direction
c
c   first relax even horizontal lines, then odd lines
c 
      do 200  ks=0,1
         iter = 0
  100    iter = iter + 1
         do 150  k=ks,kc,2
c
         kp = k + 1
         km = k - 1
c
c   fill in matrix elements for even points
c
c     point adjacent to south pole
c     using sin(phi)=sin(Phi)=-1 for boundary condition at south pole
c
      j = -jc+1
         jp = j+1
         dxds  = x(jp,k)+1.
         dxdz2 = x(j,kp)-2.*x(j,k)+x(j,km)
         f2    = x(jp,k)-1.
         f3    = 1. - 0.25*f2*f2
         f4    = (0.25*(x(jp,kp)-x(jp,km)))**2
         term2 = f1(j)*f4*(1. + 0.75*f2*f2)/(f3*f3*f3)
         u(j,k)=  dxdz2 + term2
         d(j,k)= -2.*dxds
         g(j,k)= -(dxds*dxdz2 + f1(j)*f2*f4/(f3*f3)
     2         + gamma(j,k)*sigma(j,k)*fac2)
c
c     point adjacent to north pole
c     using sin(phi)=sin(Phi)=1 for boundary condition at north pole
c
      j =  jc-1
         jm = j-1
         dxds  = 1.-x(jm,k)
         dxdz2 = x(j,kp)-2.*x(j,k)+x(j,km)
         f2    = 1.+x(jm,k)
         f3    = 1. - 0.25*f2*f2
         f4    = (0.25*(x(jm,kp)-x(jm,km)))**2
         term2 = f1(j)*f4*(1. + 0.75*f2*f2)/(f3*f3*f3)
         d(j,k)= -2.*dxds
         l(j,k)= -dxdz2 + term2
         g(j,k)= -(dxds*dxdz2 + f1(j)*f2*f4/(f3*f3)
     2         + gamma(j,k)*sigma(j,k)*fac2)
      do 105  j=-jc+3,jc-3,2
         jm = j-1
         jp = j+1
         dxds  = x(jp,k)-x(jm,k)
         dxdz2 = x(j,kp)-2.*x(j,k)+x(j,km)
         f2    = x(jp,k)+x(jm,k)
         f3    = 1. - 0.25*f2*f2
         f4    = (0.25*(x(jp,kp)+x(jm,kp)-x(jp,km)-x(jm,km)))**2
         term2 = f1(j)*f4*(1. + 0.75*f2*f2)/(f3*f3*f3)
         u(j,k)=  dxdz2 + term2
         d(j,k)= -2.*dxds
         l(j,k)= -dxdz2 + term2
         g(j,k)= -(dxds*dxdz2 + f1(j)*f2*f4/(f3*f3)
     2         + gamma(j,k)*sigma(j,k)*fac2)
  105 continue
c
c   fill in matrix elements for odd points 
c
      do 110  j=-jc+2,jc-2,2
         jm = j-1
         jp = j+1  
         xjs= x(j,k)*x(j,k)
         fj = xjs - scs(j)
         dj = 1. - xjs
         l(j,k) = -1.
         d(j,k) = fac1(j)*(2.*x(j,k)*(1.-scs(j))/(dj*dj))
         u(j,k) = 1.
         g(j,k) = -(fac1(j)*(fj/dj) + x(jp,k) - x(jm,k))
  110 continue 
  150 continue
c     print *, 'factor matrix for ks =',ks  
c
c   do matrix factorization
c
      mm = (kc + 2 - ks) / 2
      nn = jt
      inc = 1
      jump = 2*jt 
      js   = -jc + 1
      call gtdmsf(mm,nn,l(js,ks),d(js,ks),u(js,ks),inc,jump)
c
c   solve matrix system
c
c     print *, 'solve matrix for ks =',ks  
      call gtdmss(mm,nn,l(js,ks),d(js,ks),u(js,ks),g(js,ks),inc,jump)
c
c   update the solution
c
      do 170  k=ks,kc,2
      do 170  j=-jc+1,jc-1
         x(j,k) = x(j,k) + g(j,k)
  170 continue
      if(iter .lt. niter) go to 100 
  200 continue
c
c   *************** compute the residual *******************************
c
      do 205  k=0,kc
         do 205  j=-jc+1,jc-1
         rn  = rn + g(j,k)**2
  205 continue
      rn = sqrt(rn)
c
c     print out the residual norm every nup iterations
c
      if(mod(n-1,nup) .eq. 0)  then
         print 208, n, rn
      end if
  208 format('  iteration count = ',i5,' residual norm = ',1pe10.3)
c
c   ************************ update boundaries ***********************
c
c     ****** top *******
c
c     even points
c
      do 210  j=-jc+1,jc-1,2
          x(j,kcp) = x(j,kcm) + tdz*(1.-alpha)**kappa 
     2             / (kappa*alpha)
  210 continue
c
c     odd points
c
      do 220  j=-jc+2,jc-2,2
         jm = j - 1
         jp = j + 1
         x(j,kcp) = x(j,kc)
         do 215 itr=1,4
            xold     = x(j,kcp)
            xs       = xold*xold
            gkp      = fac1(j)*((xs-scs(j))/(1.-xs))
     2               + x(jp,kcp) - x(jm,kcp)
            gkpp     = fac1(j)*2.*xold*(1.-scs(j))
     2               / (1.-xs)**2
            x(j,kcp) = xold - gkp/gkpp
  215    continue
  220 continue
c
c     ****** bottom ******
c
c     even points
C
      DO 225  j=-jc+1,jc-1,2
         jp = j + 1
         jm = j - 1
         if(j .eq. -jc+1) then
            xsp = x(jp,0)*x(jp,0)
             t1 = 0.0
             t2 = (xsp-scs(jp))**2 / (1.-xsp)
         elseif(j .eq. jc-1) then
            xsm = x(jm,0)*x(jm,0)
             t1 = (xsm-scs(jm))**2 / (1.-xsm) 
             t2 = 0.0
         else
            xsp = x(jp,0)*x(jp,0)
            xsm = x(jm,0)*x(jm,0)
             t1 = (xsm-scs(jm))**2 / (1.-xsm) 
             t2 = (xsp-scs(jp))**2 / (1.-xsp)
         endif
         term3 = eps*dz*(t1+t2)/(8.*beta)
         x(j,-1) = x(j,1) - tdz/beta*x(j,0) + term3 
  225 CONTINUE
c
c     odd points
c
      do 235  j=-jc+2,jc-2,2
         jm = j - 1
         jp = j + 1
         x(j,-1) = x(j,0)
         do 230 itr=1,4
            xold     = x(j,-1)
            xs       = xold*xold
            gkm      = fac1(j)*((xs-scs(j))/(1.-xs))
     2               + x(jp,-1) - x(jm,-1)
            gkmp     = fac1(j)*2.*xold*(1.-scs(j))
     2               / (1.-xs)**2
            x(j,-1) = xold - gkm/gkmp
  230    continue
  235 continue
c
c     find gamma for general case
c
      if(gam .eq. 'g') then 
         do 245  k=0,kc
            kp = k + 1
            km = k - 1
            do 240  j=-jc+1,jc-1,2
               jp = j + 1
               jm = j - 1
               dmz        = x(j,kp) - x(j,km)
               pie        = kappa*alpha*dmz/tdz
               gamma(j,k) = pie**((kappa-1.)/kappa)
  240       continue
  245    continue
      endif
  280 continue
  290 continue
c
c   interpolate potential latitude and theta postions to actual
c   latitude
c
      call int(x, sc, convd, itm/dhr, thdot, sa, lambda)
      pmin = phimin
      pmax = phimax
      tmin = theta(0)
      tmax = theta(kc)
      call heating(itm, pmin, pmax, tmin, tmax)
c
  999 continue
c
c   plot position of parcels
c
      call traj(shr, ehr, dhr, pmin, pmax, tmin, tmax)
      go to 45
 1000 continue
c     call heating(pmin, pmax, tmin, tmax)
c
c   close workstation  
c   
      call clsgks
c
      end
      subroutine traj(shr, ehr, dhr, pmin, pmax, tmin, tmax)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy
      parameter (jl=60, kc=40)
c
c   routine to plot parcel trajectories
c
      integer shr, ehr, dhr, kinc, jinc, j, k, n 
      integer js, nhr, ns, ne
      real  x(0:8), y(0:8)
      real  pmin, pmax, tmin, tmax
      character labx*80, laby*80, title*80, ipat*5
      real  ppos, tpos, thl
      logical row
c
      common /pos/   ppos(-jl:jl,0:kc,0:8), tpos(-jl:jl,0:kc,0:8),
     2               thl(-jl:jl,0:kc)
      nhr = (ehr - shr + 1.e-5) / dhr
      ns  = shr / dhr
      ne  = ehr / dhr
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'theta (K)'
      mady = 3
      midy = 2
      madx = 6
      midx = 1
      nsdx = 0
      nsdy = 0
      write(title,60) shr/24., ehr/24.
   60 format('parcel trajectories for day ',f4.2,' to ',f4.2)
      ipat = '$$$$$' 
      call dashdc(ipat, 5, 2)
      call pfhq
      call pfdupl
      call pfsetr( 'RATIO', 0.6)  
      call setbig(pmin, pmax, tmin, tmax, 1)
   70 continue
      print 71
   71 format('Specify increment in vertical to plot trajectories',
     2        /' (0 to quit) ...>')
      read(5,*) kinc
      if(kinc .eq. 0) go to 500
      print 72
   72 format('Specify increment (even number) in latitude to',
     2         '  plot ...>')
      read(5,*) jinc
      row   = .true.
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call ptitle(title, 0.8)
      do 100  k=1,kc-1,kinc
         if(row) then
            js = -jl+1
            row   = .false.
         else
            js = -jl+1 + jinc/2
            row   = .true.
         endif
         do 100 j=js,jl-1,jinc
            do 90 n=ns,ne
               x(n) = ppos(j,k,n)
               y(n) = tpos(j,k,n)
               write(7,*) k, j, n, x(n), y(n)
   90       continue
            call curved(x(ns), y(ns), nhr)
  100 continue
      call frame
      go to 70
  500 continue
      return
      end
      subroutine int(x, sc, convd, ntim, thdot, sa, lambda)
c
c   routine to interpolate values to latitude space
c
      integer  jc, kc, kcp, jc2, jl, nwk, ntim
      parameter (kc = 40,  kcp=kc+1)
      parameter (jc = 81, jc2=jc/2)
      parameter (nwk=2*(jc+1))
      parameter (jl = 60)
      integer  j, k, jm, jp, jd2 
c
      SINGLE_OR_DOUBLE x(-jc+1:jc-1,-1:kcp), sc(-jc:jc), lambda
      SINGLE_OR_DOUBLE plm, plp, convd, thetat, thetab, pi, lz
      real  plat(-jc2:jc2,0:kc), lat, dx(-jc2:jc2), dy(0:kc)
      real  phibt(-jl:jl), d2(-jc2:jc2)      
      real  pp(-jl:jl,0:kc), th0(-jc2:jc2,0:kc)
      real  pllat(-jc2:jc2), wk(nwk), thet, yp(0:kc), ff(0:kc)
      real  thdot(-jc:jc,0:kc), thd(-jc2:jc2,0:kc), ffe(0:kc) 
c     real  tht(-jl:jl,0:kc)
      real  ppos, tpos, thl, snz, zeta, qht, thdt
      character sa
c
      common /param/ thet(0:kc), qht(-jl:jl), lat(-jl:jl), thetat,
     2               thetab, pi
      common /pos/   ppos(-jl:jl,0:kc,0:8), tpos(-jl:jl,0:kc,0:8),
     2               thl(-jl:jl,0:kc)
      save thd
c
c   compute potential latitude in degrees
c
      pllat(-jc2) = -90.
      pllat( jc2) =  90.
      do 10 j=-jc+3,jc-3,2
         pllat(j/2) = asin(sc(j))*convd
   10 continue 
      if(ntim .eq.0) then
         do 38  k=0,kc
          do 38  j=-jc,jc,2
            jd2 = j/2
            plat(jd2,k) = pllat(jd2)
            th0(jd2,k) = thet(k)
            thd(jd2,k) = thdot(j,k)
   38    continue
         call spx(plat, pllat, jc2, kc, lat, jl, wk, nwk, dx,  
     2            ppos(-jl,0,0))
         call spx(th0, pllat, jc2, kc, lat, jl, wk, nwk, dx,  
     2            tpos(-jl,0,0))
         call spx(thd,  pllat, jc2, kc, lat, jl, wk, nwk, dx, thl)
      else
      do 20 k=0,kc
         plat(-jc2,k) = -90.
         plat( jc2,k) =  90.
         do 15  j=-jc+3,jc-3,2
            jm = j - 1
            jp = j + 1
            plm = asin(x(jm,k))
            plp = asin(x(jp,k))
            plat(j/2,k) = 0.5*(plp+plm)*convd
   15    continue
   20 continue
c
c   fit spline phil_pos(phi_o) at each theta_o
c 
         call spx(plat, pllat, jc2, kc, lat, jl, wk, nwk, dx, pp) 
c
c   fit spline phil_pos(theta_o) at each phi_o
c   and evaluate at theta(ntim-1)
c
         call spy(pp, thet, jl, kc, tpos(-jl,0,ntim-1), kc, wk, nwk,
     2         dy, ppos(-jl,0,ntim), yp, ff, ffe)
c 
c   compute heating field in little phi space
c
         call heat(thd, plat, pllat, jc2, kc, 
     2             lat, phibt, jl, wk, nwk, dx, d2, thl)
c
c   compute theta position
c
         do 70 k=0,kc
            do 60 j=-jl,jl
               zeta = (tpos(j,k,ntim-1) - thetab)/(thetat-thetab)
               lz = lambda * zeta
               if(lz .gt. 1) then
                  thdt = 0.0
               else
                  snz  = sin(pi*lz)
                  if(sa .eq. 'b') then
                     thdt = qht(j)*snz
                  elseif(sa .eq. 'i') then
                     thdt = qht(j)*snz*snz
                  endif
               endif
               tpos(j,k,ntim) = tpos(j,k,ntim-1) + thdt
   60       continue
   70    continue
      endif
      return
      end
      subroutine heat(thd, phil, phib, jj, kk, 
     2                lat, phibt, je, wk, nwk, d1, d2, thl)
c
c   routine to spline fit a field (f) as a function of independent
c   variable (x) and then evaluate the field at points xe and
c   put answer in fe.
c
      real  thd(-jj:jj,0:kk), phil(-jj:jj,0:kk), phib(-jj:jj)
      real  d1(-jj:jj), d2(-jj:jj), vc(2), switch
      real  lat(-je:je), phibt(-je:je), thl(-je:je,0:kk), wk(nwk)
      integer ic(2), n, nwk, ne, ierr, k, je, kk, jj
      logical skip
c
      ic(1) = -1
      vc(1) = 0.0
      ic(2) = -1
      vc(2) = 0.0
      switch = -1
      skip = .true.
      n =  2*jj + 1
      ne = 2*je + 1
      do 20  k=0,kk
c   fit spline phib(phil)
         call pchic( ic, vc, switch, n, phil(-jj,k), phib,
     2               d1, 1, wk, nwk, ierr)
c   fit spline thd(phib)
         call pchic( ic, vc, switch, n, phib, thd(-jj,k),
     2               d2, 1, wk, nwk, ierr)
c   evaluate phib(phil) at lat points >> phibt.
         call pchfe( n, phil(-jj,k), phib, d1, 1, skip,
     2               ne, lat, phibt, ierr) 
c   evaluate thd(phibt) at lat points >> thl.
         call pchfe( n, phib, thd(-jj,k), d2, 1, skip,
     2               ne, phibt, thl(-je,k), ierr) 
   20 continue
      return
      end
      subroutine spx(f, x, jj, kk, xe, je, wk, nwk, d, fe)
c
c   routine to spline fit a field (f) as a function of independent
c   variable (x) and then evaluate the field at points xe and
c   put answer in fe.
c
      real  f(-jj:jj,0:kk), x(-jj:jj)
      real  d(-jj:jj), vc(2), switch
      real  xe(-je:je), fe(-je:je,0:kk), wk(nwk)
      integer ic(2), n, nwk, ne, ierr, k, je, kk, jj
      logical skip
c
      ic(1) = -1
      vc(1) = 0.0
      ic(2) = -1
      vc(2) = 0.0
      switch = -1
      skip = .true.
      n =  2*jj + 1
      ne = 2*je + 1
      do 20  k=0,kk
         call pchic( ic, vc, switch, n, x, f(-jj,k),
     2               d, 1, wk, nwk, ierr)
         call pchfe( n, x, f(-jj,k), d, 1, skip,
     2               ne, xe, fe(-je,k), ierr) 
   20 continue
      return
      end
      subroutine spy(f, y, jj, kk, ye, ke, wk, nwk, d, fe, 
     2               yp, ff, ffe)
c
c   routine to spline fit a field (f) as a function of independent
c   variable (y) and then evaluate the field at points ye and
c   put answer in fe.
c
      real  f(-jj:jj,0:kk), y(0:kk), ff(0:kk), yp(0:ke), ffe(0:ke)
      real  d(0:kk), vc(2), switch
      real  ye(-jj:jj,0:ke), fe(-jj:jj,0:ke), wk(nwk)
      integer ic(2), n, nwk, ne, ierr, j, jj, kk, ke, k
      logical skip
c
      ic(1) = -1
      vc(1) = 0.0
      ic(2) = -1
      vc(2) = 0.0
      switch = -1
      skip = .true.
      n = kk + 1
      ne = ke + 1
      do 20  j=-jj,jj
         do 15  k=0,ke
           ff(k) =  f(j,k)
           yp(k) = ye(j,k)
   15    continue
         call pchic( ic, vc, switch, n, y, ff,
     2               d, 1, wk, nwk, ierr)
         call pchfe( n, y, ff, d, 1, skip,
     2               ne, yp, ffe, ierr) 
         do 17  k=0,ke
           fe(j,k) =  ffe(k)
c          if(j.eq.0) print *, j, k, y(k), ff(k), yp(k), ffe(k)
   17    continue
   20 continue
      return
      end
      subroutine heating(itm, pmin, pmax, tmin, tmax)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy, jtot, itm
      parameter (jl=60, kc=40, jtot=2*jl+1)
c
c   routine to plot heating
c
      real  pmin, pmax, tmin, tmax
      character labx*80, laby*80, title*80
      real  pposl, tposl, thl
c
      common /pos/   pposl(-jl:jl,0:kc,0:8), tposl(-jl:jl,0:kc,0:8),
     2               thl(-jl:jl,0:kc)
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'theta (K)'
      write(title, 5) itm/24.
    5 format(' heating (K/day) at ',f4.2,' days')
      mady = 3
      midy = 2
      madx = 6
      midx = 1
      nsdx = 0
      nsdy = 0
      call pfhq
      call pfdupl
      call pfsetr( 'RATIO', 0.6)  
      call setbig(pmin, pmax, tmin, tmax, 1)
      call pcseti('QU - quality flag for labels', 0)
c     call cpseti('NCL - number of contour levels ', 5)
      call cpseti('CLS - contour selection flag ', 5)
      call cpsetr('CIS - contour interval specifier',1.)
      call cpseti('LIS - line interval specifier',1)
      call cpsetr('PC1 - penalty scheme constant 1', 20.)
      call cpsetr('PC3 - penalty scheme constant 3', 200.)
      call cpseti('LLP - line label positioning', 3)
      call cpsetr('ILY - info label y-coor.', 1.2)
      call cpseti('T2D - tension on 2d smoother', 4)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call ptitle(title, 0.8)
      call cpcnrc(thl, jtot, jtot, kc+1, 0.0, 0.0, 1.,
     2            -1, -1, -924)
      call frame
      return
      end
