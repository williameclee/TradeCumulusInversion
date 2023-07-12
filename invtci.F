      program invert
c
c
c     This program uses a newton iteration to solve the invertibility
c     equation for the trade cumulus inversion problem.
c
      integer  jc, jt, kc, kcp, kcm, jl
      parameter (kc = 60,  kcp=kc+1, kcm=kc-1)
      parameter (jc = 81, jt=2*jc-1, jl=60)
      integer  jm, jp, j, k, iter, nmax, n, ks, kp, km, itr
      integer  nn, inc, jump, js, mm, nup, niter, method
      real  x(-jc+1:jc-1,-1:kcp), l(-jc+1:jc-1,0:kc)
      real  d(-jc+1:jc-1,0:kc), u(-jc+1:jc-1,0:kc)
      real  g(-jc+1:jc-1,0:kc), phid(-jc:jc)
      real  sc(-jc:jc), scs(-jc:jc)
      real  fac1(-jc:jc), f1(-jc:jc)
      real  z(0:kc), theta(0:kc), pres(-jc:jc,0:kc)
      real  gamma(-jc:jc,0:kc), sigma(-jc:jc,0:kc)
      real  cs, eps, ae, pt, pb, sigma0, avepb 
      real  r, gr, omega, rom, cp, kappa, avep, dif
      real  pi, pid2, sqpi, dz, tdz, dt, alpha, beta
      real  convd, convr, sum
      real  xjs, fj, dj, dsc, rn, refp(0:kc)
      real  dxds, dxdz2, f2, f3, f4, term2
      real  xold, fac2, gkp, gkpp, xs, xsp, xsm, term3
      real  gkm, gkmp, dmz, pie, phimax, phimin
      real  t1, t2, ftmp(0:kc), xtmp(0:kc)
      real              theta10, thetac, theta1(-jc:jc), theta2(-jc:jc)
      real              tm, tp, dfdt, fphi, fthet, tinc, fm, fp
      real              thetab, thetat, dp, phis, phin, phic
      real              dtheta, fexp, lat(-jl:jl), dl, ptop
      common /coord/ theta, lat
      common /param/ kappa, alpha, tdz, pb, pt, omega, ae, convr,
     2               cs, cp, gr 
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
      pid2  = pi / 2.d0
      convd = 180.d0 / pi
      convr = pi / 180.d0
      sqpi  = sqrt(pi) 
c
c   read in theta at model top (thetat) and bottom (thetab)
c
      print *, 'Specify theta (kelvin) at model top      ...>'
      read(5,*) thetat
      print *, 'Specify theta (kelvin) at model bottom   ...>'
      read(5,*) thetab
      print *, 'Specify pressure (mb) at model bottom    ...>'
      read(5,*) pb
      print *, 'Specify inversion pressure change (mb)   ...>'
      read(5,*) dp
      print *, 'Specify sigma0 (kPa/K)                   ...>'
      read(5,*) sigma0
      print *, 'Specify mean theta (kelvin) in inversion ...>'
      read(5,*) theta10
      print *, 'Specify dtheta (kelvin) of inversion     ...>'
      read(5,*) dtheta
      print *, 'Specify thetac (kelvin) of inversion     ...>'
      read(5,*) thetac
      print *, 'Specify southern extent of inversion' 
      print *, ' latitude (in degrees)                   ...>'
      read(5,*) phis
      print *, 'Specify northern extent of inversion' 
      print *, ' latitude (in degrees)                   ...>'
      read(5,*) phin
      print *, 'Specify center of inversion' 
      print *, ' latitude (in degrees)                   ...>'
      read(5,*) phic
c
c   compute pt in units of Pa 
c
      dp = dp * 100.
      pb = pb * 100.
      sigma0 = sigma0 * 1000.
      write(7,*) sigma0
      pt = pb - sigma0*(thetat-thetab) + dp
c
c   compute other parameters
c
      dt     = (thetat - thetab) / kc
      alpha  = (pb-pt) / pb
      beta   = thetab / (thetat-thetab)
      print *, ' thetab (K) = ', thetab
      print *, ' thetat (K) = ', thetat
      print *, ' pb     (mb) = ', pb/100.
      print *, ' pt     (mb) = ', pt/100.
      print *, ' sigma0 (Pa/K) = ', sigma0
c
      cs    = alpha*r*(thetat-thetab)
      eps   = 4.d0 * rom*rom / cs
      print 2, eps, sqrt(cs)
    2 format(' epsilon = ',1p,e10.3,' c = ',e10.3)
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
   20 continue
c
c     ****************** INITIALIZE VARIABLES *******************
c
c   initialize theta1 and theta2
c
      do 25  j=-jc,jc
         if(phid(j).gt.phin .or. phid(j).lt.phis) then
            theta1(j) = theta10
         elseif(phid(j).le.phin .and. phid(j).gt.phic) then
            fphi      = (phin-phid(j)) / (phin-phic)
            theta1(j) = theta10 + dtheta * fexp(fphi) 
         elseif(phid(j) .eq. phic) then
            theta1(j) = theta10 + dtheta 
         elseif(phid(j).lt.phic .and. phid(j).ge.phis) then
            fphi      = (phid(j)-phis) / (phic-phis)
            theta1(j) = theta10 + dtheta * fexp(fphi) 
         endif           
         theta2(j) = theta1(j) + thetac
c        write(7,*) j, phid(j), theta1(j), theta2(j)
   25 continue
c           
c   initialize pressure, gamma and x (i.e., m) at even points
c
      do 30  k=0,kc
         do 28  j=-jc,jc
            if(theta(k) .ge. theta2(j)) then
               pres(j,k) = pb - sigma0*(theta(k)-thetab)
     2                        + dp
            elseif(theta(k).lt.theta2(j) .and. 
     2             theta(k).ge.theta1(j)) then
               fthet = (theta(k)-theta1(j)) / (theta2(j)-theta1(j))
               pres(j,k) = pb - sigma0*(theta(k)-thetab)  
     2                        + dp*fexp(fthet)
            elseif(theta(k) .lt. theta1(j)) then
               pres(j,k) = pb - sigma0*(theta(k)-thetab)  
            endif
            pres(j,k)  = pres(j,k)/pb 
            gamma(j,k) = pres(j,k)**(kappa-1.) 
   28    continue
         sum = 0.0
         do 29  j=-jc+1,jc-1,2
            sum = sum + pres(j,k)*(sc(j+1)-sc(j-1))
   29    continue
         refp(k) = 0.5 * sum 
c        write(7,*) k, refp(k)*pb
   30 continue
c
c     print out the residual norm every nup iterations
c
c   compute initial m field at even points
c
      method = 1
      do 40 j=-jc+1,jc-1,2
         x(j,0) = cp*thetab/cs
         do 35  k=0,kc
            ftmp(k) = pres(j,k)**kappa/(kappa*alpha)  
   35    continue
         call intdde( method, kc, dz, x(j,0), ftmp, xtmp)
         do 37  k=0,kc
            x(j,k) = xtmp(k)
   37    continue
         x(j,-1)  = x(j,1) - tdz*x(j,0)/beta
         x(j,kcp) = x(j,kcm) + tdz*(1.d0-alpha)**kappa 
     2                       / (kappa*alpha)
c        if(j .eq. 0) then
c           do 39  k=-1,kcp
c              write(7,*) k, x(0,k)
c  39       continue
c        endif
   40 continue
c
c     print out the residual norm every nup iterations
c
c   initialize x (i.e., sin(phil)) at odd points
c
      do 48  k=-1,kcp
      do 48  j=-jc+2,jc-2,2
         x(j,k) = sc(j)
   48 continue
c
c     ******************* COMPUTE SIGMA ***********************
c
      tinc = 0.1
      DO 70 k=0,kc
         DO 60  j=-jc,jc
            if(theta(k).ge.theta2(j) .or. theta(k).le.theta1(j)) then
               sigma(j,k) = 1.
            elseif(theta(k).lt.theta2(j) .and. 
     2             theta(k).ge.theta1(j)) then
               tp   = theta(k) + tinc
               tm   = theta(k) - tinc
               fp   = (tp-theta1(j)) / (theta2(j)-theta1(j))
               fm   = (tm-theta1(j)) / (theta2(j)-theta1(j))
               dfdt = (fp - fm) / (2.*tinc)
               sigma(j,k) = (sigma0 - dp*dfdt)/sigma0      
            endif
   60    continue
   70 continue
      print *, 'Specify phimin (degrees) for plot',
     2         '(e.g., phimin = -30, phimax =0) ...>'
      read(5,*) phimin
      print *, 'Specify phimax (degrees) for plot...>'
      read(5,*) phimax
      print *, 'Phimin and Phimax are: ', phimin, phimax
      print *, 'Specify ptop (kPa) for pressure plot...>'
      read(5,*) ptop   
      call plot(x, sigma, theta, z, phimin, phimax, scs, sc, 
     2          ptop, 'initial')
c
c   interpolate potential latitude and theta postions to actual
c   latitude
c
      dl = (phimax - phimin) / (2*jl)
      do 54  j=-jl,jl
         lat(j) = phimin + (j+jl)*dl
   54 continue
      call int(x, sc, convd, 0)
c     go to 999
C
C     ******************* BEGIN ITERATION *********************
c
c     specify nmax (the number of relaxation sweeps to be performed)
c
      print *,'Specify number of sweeps  ...>'
      read  *, nmax
      print *,'Number of sweeps  to be done: ',nmax
      print *,'Specify how frequent to print residual ...>'
      read  *, nup
      print *,'Frequency to print residual at: ',nup
c   do "niter" newton iteration at every sweep
      niter = 1
      n  = 0
   80 n  = n + 1
      if (n .gt. nmax) go to 300
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
C
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
  225 continue
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
      do 245  k=0,kc
         kp = k + 1
         km = k - 1
         sum = 0.0
         do 240  j=-jc+1,jc-1,2
            dmz        = x(j,kp) - x(j,km)
            pie        = kappa*alpha*dmz/tdz
            pres(j,k)  = pie**(1./kappa)
            if(j .eq. -jc+1) then
               sum = sum + pres(j,k)*(x(j+1,k)+1.)
            else if(j .eq. jc-1) then
               sum = sum + pres(j,k)*(1.-x(j-1,k))
            else
               sum = sum + pres(j,k)*(x(j+1,k)-x(j-1,k))
            endif
c           sum = sum + pres(j,k)*(sc(j+1)-sc(j-1))
  240    continue
         avep = 0.5 * sum 
         dif  = avep - refp(k)
         if(k .eq. 0) avepb = avep
         sum = 0.0
         do 242  j=-jc+1,jc-1,2
            pres(j,k)  = pres(j,k) - dif
            if(j .eq. -jc+1) then
               sum = sum + pres(j,k)*(x(j+1,k)+1.)
            else if(j .eq. jc-1) then
               sum = sum + pres(j,k)*(1.-x(j-1,k))
            else
               sum = sum + pres(j,k)*(x(j+1,k)-x(j-1,k))
            endif
c           sum = sum + pres(j,k)*(sc(j+1)-sc(j-1))
            gamma(j,k) = pres(j,k)**(kappa-1.)
  242    continue
         avep = 0.5 * sum 
         dif  = avep - refp(k)
c        write(7,*) 'iteration = ', n, avep, refp(k), dif
  245 continue
c
c   compute a corrected m field 
c
c     go to 261
      do 260 j=-jc+1,jc-1,2
         do 255  k=0,kc
            ftmp(k) = pres(j,k)**kappa/(kappa*alpha)  
  255    continue
         call intdde( method, kc, dz, x(j,0), ftmp, xtmp)
         do 257  k=0,kc
            x(j,k) = xtmp(k)
  257    continue
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
c        term3 = 0.0
         x(j,-1)  = x(j,1) - tdz/beta*x(j,0) + term3 
         x(j,kcp) = x(j,kcm) + tdz*(1.d0-alpha)**kappa 
     2                       / (kappa*alpha)
c        if(j .eq. 0) then
c           do 259  k=-1,kcp
c              write(7,*) k, x(0,k)
c 259       continue
c        endif
  260 continue
  261 continue
c
c     print out the residual norm every nup iterations
c
      if(mod(n-1,nup) .eq. 0)  then
         print 290, n, rn, avepb
      end if
  290 format(' iteration = ',i5,' res. norm = ',1pe10.3,
     2       ' average pbot = ',e10.3)
      go to 80
  300 continue
c
c   plot out fields
c
c   interpolate potential latitude and theta postions to actual
c   latitude
c
      call int(x, sc, convd, 1)
c
c   plot position of parcels
c
c     call traj(phimin, phimax, thetab, thetat)
c 
      call plot(x, sigma, theta, z, phimin, phimax, scs, sc, 
     2          ptop, 'final')
 999  continue
C
c   close workstation  
c   
      call clsgks
c
      end
      subroutine int(x, sc, convd, ntim)
c
c   routine to interpolate values to latitude space
c
      integer  jc, kc, kcp, jc2, jl, nwk, ntim
      parameter (kc = 60,  kcp=kc+1)
      parameter (jc = 81, jc2=jc/2, jl=60)
      parameter (nwk=2*(jc+1))
      integer  j, k, jm, jp, jd2
c
      real x(-jc+1:jc-1,-1:kcp), sc(-jc:jc)
      real plm, plp, convd, theta(0:kc)
      real  plat(-jc2:jc2,0:kc), lat(-jl:jl), dx(-jc2:jc2)
      real  pllat(-jc2:jc2), wk(nwk)
c     real  tht(-jl:jl,0:kc)
      real  ppos, tpos
c
      common /coord/ theta, lat
      common /pos/   ppos(-jl:jl,0:kc,0:1), tpos(-jl:jl,0:kc,0:1)
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
   38    continue
         call spx(plat, pllat, jc2, kc, lat, jl, wk, nwk, dx,
     2            ppos(-jl,0,0))
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
   15       continue
   20    continue
c
c   fit spline phil_pos(phi_o) at each theta_o
c
         call spx(plat, pllat, jc2, kc, lat, jl, wk, nwk, dx, 
     2            ppos(-jl,0,ntim) )
      endif
c
c   compute theta position
c
      do 70 k=0,kc
      do 60 j=-jl,jl
         tpos(j,k,ntim) = theta(k)                 
   60 continue
   70 continue
      return
      end
      subroutine traj(phimin, phimax, tmin, tmax)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy
      parameter (jl=60, kc=60)
c
c   routine to plot parcel trajectories
c
      real phimin, phimax
      integer kinc, jinc, j, k, n, js
      integer kumx, kumy
      real  x(0:1), y(0:1)
      real  pmin, pmax, tmin, tmax
      character labx*80, laby*80, title*80, ipat*5
      real  ppos, tpos
      logical row
c
      common /pos/   ppos(-jl:jl,0:kc,0:1), tpos(-jl:jl,0:kc,0:1)
      pmin = phimin
      pmax = phimax
      print *, 'phimin = ',pmin,' phimax = ',pmax
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
      write(title,60) 
   60 format('parcel trajectories ')
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
      call vvsetr('ZMX', 1.e5)
      do 100  k=1,kc-1,kinc
         if(row) then
            js = -jl+1
            row   = .false.
         else
            js = -jl+1 + jinc/2
            row   = .true.
         endif
         do 100 j=js,jl-1,jinc
            do 90 n=0,1
               x(n) = ppos(j,k,n)
               y(n) = tpos(j,k,n)
c              write(7,*) k, j, n, x(n), y(n)
   90       continue
c        call curved(x, y, 2)
         call drwvec(kumx(x(0)), kumy(y(0)), kumx(x(1)), kumy(y(1)),
     2               '', 0)
  100 continue
      call frame
      go to 70
  500 continue
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
