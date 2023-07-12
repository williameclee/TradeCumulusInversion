      subroutine plot(x, sigma, theta, z, phimin, phimax, 
     2                scs, sc, ptop, time)
      integer  jc, kc, kcp, jc2, jtot, jl, nwk, nwkz, ne
      parameter (kc = 60,  kcp=kc+1)
      parameter (jc = 81, jc2=jc/2, jtot=2*jc2+1)
      parameter (jl = 60, nwk=2*jtot, nwkz=2*kcp, ne=100)
      integer  jmin, jmax, jplt, j, jcm, jcp, k
      integer  jm, jp, km, kp, jpt
      character time*(*)
c
      real x(-jc+1:jc-1,-1:kcp), theta(0:kc), z(0:kc)
      real scs(-jc:jc), sc(-jc:jc), sigma(-jc:jc,0:kc)
      real phimin, phimax, scmin, scmax, pie
      real kappa, alpha, tdz, pb, pt, omega, ae, convr
      real plm, plp, sls, up, um, dx
      real cs, cp, gr
      real  zetk, sj, cj, cjp, cjm    
      real  plat(-jc2:jc2,0:kc), dp(-jc2:jc2), lat(-jl:jl)
      real  u(-jc2:jc2,0:kc), p(-jc2:jc2,0:kc)
      real  m(-jc2:jc2,0:kc)
      real  sig(-jc2:jc2,0:kc), pv(-jc2:jc2,0:kc)
      real  dpv(-jc2:jc2,0:kc), dpvl(-jl:jl,0:kc)
      real  ul(-jl:jl,0:kc), pl(-jl:jl,0:kc)
      real  ml(-jl:jl,0:kc)
      real  sigl(-jl:jl,0:kc), pvl(-jl:jl,0:kc)
      real  thet2(-jl:jl,0:kc), thetp(-jl:jl,0:ne)
      real  ptop, pbot, pinc(0:kc), finc(0:kc), d(0:kc)
      real  pe(0:ne), fe(0:ne)
      real  wk(nwk), wkz(nwkz), dl, dt
      real  tdt, dm, temp, hgt
      real pmin, pmax, tb, tt, lats, latmin, latdif
c
      common /param/ kappa, alpha, tdz, pb, pt, omega, ae, convr,
     2               cs, cp, gr

c
      cp = 1004.
      gr = 9.81
c
c   compute u(m/s) and potential latitude in degrees
c
      do 20 k=0,kc
         u(-jc2,k) = 0.
         u( jc2,k) = 0.
         plat(-jc2,k) = -90.
         plat( jc2,k) =  90.
         do 15  j=-jc+3,jc-3,2
            jm = j - 1
            jp = j + 1
            plm = asin(x(jm,k))
            sls= x(jm,k)*x(jm,k)
            up = omega*ae*(sls-scs(jm))/cos(plm)
            plp = asin(x(jp,k))
            sls= x(jp,k)*x(jp,k)
            um = omega*ae*(sls-scs(jp))/cos(plp)
            u(j/2,k) = (up + um) / 2.
            plat(j/2,k) = 0.5*(plp+plm)/convr
   15    continue
   20 continue
   21 format(f6.1,i5,1x,f7.2,',',1x,f7.2)
c
c   compute pressure
c
      jcm = (-jc+1)/2
      jcp = ( jc-1)/2
      do 30 k=0,kc
         kp = k + 1
         km = k - 1
         do 25 j=-jc+1,jc-1,2
            dx  = x(j,kp) - x(j,km)
            pie = kappa*alpha*dx/tdz
            p(j/2,k) = pb*pie**(1./kappa)/1000. 
            m(j/2,k) = (x(j,k) - 0.5*u(j/2,k)*u(j/2,k)) * cs
c           if(j.eq.0) write(7,*) k, theta(k), p(j/2,k)
   25    continue
         p(-jc2,k) = p(jcm,k)
         p( jc2,k) = p(jcp,k)
   30 continue
c
c   compute non-dimesional pv
c
      do 38  k=0,kc
      do 38  j=-jc+1,jc-1,2
         sig(j/2,k) = sigma(j,k)
         pv(j/2,k) = sc(j) / sig(j/2,k)
   38 continue
c
c   compute the derivative of pv field
c
      do 40 k=0,kc
         do 39 j=-jc2+1,jc2-1
            jp = j + 1
            jm = j - 1
            dpv(j,k) = pv(jp,k) - pv(jm,k)
   39    continue
         dpv(-jc2,k) = dpv(-jc2+1,k)
         dpv( jc2,k) = dpv( jc2-1,k)
   40 continue
c
c   interpolate fields to actual latitude
c
      dl = (phimax - phimin) / (2*jl)
      do 50  j=-jl,jl
         lat(j) = phimin + (j+jl)*dl
   50 continue
c
      call xtranf(u, ul, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
      call xtranf(p, pl, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
      call xtranf(m, ml, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
      call xtranf(pv, pvl, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
      call xtranf(sig, sigl, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
      call xtranf(dpv, dpvl, jc2, kc, plat, lat, jl, wk, nwk, dp, jtot)
c
c   save theta(pres) data at a specified latitude
c
   44 continue
      write(6,45) phimax, phimin
   45 format('Specify a latitude between ',f5.1,' and ',f5.1,' to save',
     2       /,'theta(pres) data for (outside this range to quit)...>')
      read(5,*) lats
      write(7,*) lats
      print *, 'saving data at ',lats
      if(lats.le.phimax .and. lats.ge.phimin) then 
         latmin = 1.e10

         do 46 j=-jl,jl
            latdif = abs(lat(j) - lats)
            if(latdif .lt. latmin) then
               latmin = latdif
               jpt    = j
            endif
   46    continue

c        write(7,47) lat(jpt)
   47    format(/,'theta(K) vs pressure(mb) at lat =',f7.2)
c        cj = cos(lat(jpt)*convr)
c        cjp = cos(lat(jpt+1)*convr)
c        cjm = cos(lat(jpt-1)*convr)
c        sj  = sin(lat(jpt)*convr)

         tdt = 2.*(theta(1)-theta(0))

         do 48  k=0,kc  
c           zetk = -(ul(jpt+1,k)*cjp - ul(jpt-1,k)*cjm)/
c    2                          (ae*cj*(lat(jpt+1)-lat(jpt-1)))
c           write(7,49) k, theta(k), pl(jpt,k)*10., zetk, pvl(jpt,k)
            if(k .eq. 0) then
               dm = -3.*ml(jpt,k) + 4.*ml(jpt,k+1) - ml(jpt,k+2)
            elseif(k .eq. kc) then
               dm =  3.*ml(jpt,k) - 4.*ml(jpt,k-1) + ml(jpt,k-2)
            else
               dm = ml(jpt,k+1) - ml(jpt,k-1)
            endif
            pie = dm / tdt 
            temp = pie * theta(k)/cp
            hgt  = (ml(jpt,k) - theta(k)*pie) / gr
            write(7,49) k, theta(k), pl(jpt,k)*10., temp, hgt                        

   48    continue

   49    format(i3,4e14.6)
         go to 44
      endif
c
c   plot out fields in latitude space
c
      pmin = phimin
      pmax = phimax
      tb   = theta(0)
      tt   = theta(kc)
c
c   plot zonal winds (u) and pressure (p) in latitude space
c
      call pltlup(ul, pl, jl, kc, pmin, pmax, tb, tt, time)
c
c   plot sigma in latitude space
c
      call pltls(sigl, jl, kc, pmin, pmax, tb, tt, time)
c
c   plot potentail vorticity and shade area of neg. pv gradient
c
      call pltlpv(pvl, dpvl, jl, kc, pmin, pmax, tb, tt, time)
c
c   plot fields in latitude, pressure space
c
      do 55  k=0,kc
      do 55  j=-jl,jl
         thet2(j,k) = theta(k)
   55 continue 
c     pbot = 90.
      pbot = pb/1000. 
      call ptranf(thet2, thetp, jl, kc, pl, pinc, finc, d,
     2            pe, fe, ne, wkz, nwkz, pbot, ptop)
      call pltlp(thetp, jl, ne, pmin, pmax, pbot, ptop, time) 
c
c     call ptranf(u, up, jl, kc, pl, pinc, finc, d,
c    2            pe, fe, ne, wk, nwkz, pbot, ptop)
c
c   compute limits for plot: jmin and jmax
c
      scmin = sin(phimin*convr)
      scmax = sin(phimax*convr)
      call limits(scmin, scmax, sc, jc, convr, jmin, jmax, jplt)
c
c   draw plot
c
      pmin = scmin
      pmax = scmax
c     call pltpl(u,p,jc2,jtot,kc,pmin,pmax,tb,tt,time,jmin,jplt)
      return
      end
      subroutine xtranf(fpl, fl, jc2, kc, plat, lat, jl,
     2           wk, nwk, dp, jtot)
c
c   routine to transform fields from potential latitude to actual   
c   latitude
c   fpl - field in potential latitude space
c   fl  - field in latitude sapce
c
      real  vc(2), switch, dp(-jc2:jc2), wk(nwk)
      real  fl(-jl:jl,0:kc), lat(-jl:jl)
      real  fpl(-jc2:jc2,0:kc), plat(-jc2:jc2,0:kc)
      integer ic(2), n, nwk, ne, ierr, k, jc2, kc, jl, jtot
      logical skip
c
      ic(1) = -1
      vc(1) = 0.0
      ic(2) = -1
      vc(2) = 0.0
      switch = -1
      skip = .true.
      n = jtot
      ne = 2*jl + 1
      do 20  k=0,kc
         call pchic( ic, vc, switch, n, plat(-jc2,k), fpl(-jc2,k),
     2               dp, 1, wk, nwk, ierr)
         call pchfe( n, plat(-jc2,k), fpl(-jc2,k), dp, 1, skip,
     2               ne, lat, fl(-jl,k), ierr) 
   20 continue
      return
      end
      subroutine ptranf(ft, fp, jl, kc, pres, pinc, finc, d,
     2                  pe, fe, ne, wk, nwk, pb, pt)
c
c   routine to transform field from theta to pressure space
c   fp - field in pressure space          
c   ft - field in theta sapce
c
      real  vc(2), switch, wk(nwk), pb, pt
      real  pres(-jl:jl,0:kc), ft(-jl:jl,0:kc) 
      real  pinc(0:kc), finc(0:kc), d(0:kc)
      real  fp(-jl:jl,0:ne), pe(0:ne), fe(0:ne)
      integer ic(2), n, nwk, ne, ierr, k, j, jl, kc, kk, l
      logical skip
c
      ic(1) = 0
      vc(1) = 0.0
      ic(2) = 0
      vc(2) = 0.0
      switch = -1.0
      skip = .true.
      n = kc + 1
      do 10 l=0,ne
         pe(l) = pb - l*(pb-pt)/ne
   10 continue
      do 50  j=-jl,jl   
         do 20  k=0,kc
            kk = kc - k
            pinc(k) = pres(j,kk)
            finc(k) = ft(j,kk)
            if(k.gt.0 .and. pinc(k).le.pinc(k-1)) then     
               pinc(k) = pinc(k-1) + 1.e-5
               write(7,*) 'had to fix pinc for ',j,k
            endif
   20    continue
         call pchic(ic,vc,switch,n,pinc,finc,d,1,wk,nwk,ierr)
         call pchfe(n,pinc,finc,d,1,skip,ne+1,pe,fe,ierr)
         do 30  l=0,ne
            fp(j,l) = fe(l)
   30    continue
   50 continue
      return
      end
