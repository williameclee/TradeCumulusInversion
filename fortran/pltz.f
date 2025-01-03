      program pltz    
c
c   plots vertical profiles of various fields generated by invtci
c   code
c
      integer    kc    
      parameter (kc = 60)    
      real   ti(0:kc), pi(0:kc), thetai(0:kc), zetri(0:kc), pvi(0:kc)
      real   tf(0:kc), pf(0:kc), thetaf(0:kc), zetrf(0:kc), pvf(0:kc)
      real   sigmai(0:kc), sigmaf(0:kc)
      real   r, cp, kappa, sigma0, lat
      character  title*(80), fname*11, labx*20, pltloc*2                  
c
c   *****************************************************************
c
c   call plotting routines
c
      r = 287.
      cp = 1004.
      omega = 7.292e-5
      kappa = r / cp
c
c   open gks window
c
      call opngks
c
c   turn off clipping
c
      call gsclip ( 0 )
      print *, 'Specify input file name ...>'
      read(5,'(a11)') fname
      open(2, file=fname, form='formatted', status='unknown')
      read(2,*) pmin   
      read(2,*) sigma0
   10 continue
      read(2,*,end=999) lat
      do 30 k=0,kc
         read(2,*) kk, thetai(k), pi(k), zetri(k), pvi(k)   
   30 continue
      do 50  k=0,kc
         ti(k) = (thetai(k) / (1000./pi(k))**kappa) - 273.15
         zetri(k) = zetri(k) / (2.*omega) * 1.e5
         if(k .eq. 0) then
            dp = -3.*pi(k)     + 4.*pi(k+1)     - pi(k+2)
            dt = -3.*thetai(k) + 4.*thetai(k+1) - thetai(k+2)
         elseif(k .eq. kc) then
            dp =  3.*pi(k)     - 4.*pi(k-1)     + pi(k-2)
            dt =  3.*thetai(k) - 4.*thetai(k-1) + thetai(k-2)
         else
            dp = pi(k+1) - pi(k-1)
            dt = thetai(k+1) - thetai(k-1)
         endif
         sigmai(k) = -dp/dt / sigma0 * 100.
   50 continue
      read(2,*) lat
      do 130 k=0,kc
         read(2,*) kk, thetaf(k), pf(k), zetrf(k), pvf(k)   
  130 continue
      do 150  k=0,kc
         tf(k) = (thetaf(k) / (1000./pf(k))**kappa) - 273.15
         zetrf(k) = zetrf(k) / (2.*omega) * 1.e5
         print *, k, zetrf(k), zetri(k)
         if(k .eq. 0) then
            dp = -3.*pf(k)     + 4.*pf(k+1)     - pf(k+2)
            dt = -3.*thetaf(k) + 4.*thetaf(k+1) - thetaf(k+2)
         elseif(k .eq. kc) then
            dp =  3.*pf(k)     - 4.*pf(k-1)     + pf(k-2)
            dt =  3.*thetaf(k) - 4.*thetaf(k-1) + thetaf(k-2)
         else
            dp = pf(k+1) - pf(k-1)
            dt = thetaf(k+1) - thetaf(k-1)
         endif
         sigmaf(k) = -dp/dt / sigma0 * 100.
  150 continue
c
c   plot fields              
c
c     write(title,200) lat
c 200 format('temperature (C) at ',f4.0,
c    2       ' initial(solid); final(dashed)')
c     call plts(ti, tf, pi, pf, kc+1, title, 0., 30., pmin)
c
      labx = 'temperature (K)'
      pltloc = 'ul'
      write(title,210) lat
  210 format('theta (K) at ',f4.0)
      call plts(thetai, thetaf, pi, pf, kc+1, title, 300., 315., 
     2          pmin, labx, pltloc)
c
      labx = 'sigma/sigma0'
      pltloc = 'll'
      write(title,220) lat
  220 format('sigma/sigma0 at ',f4.0)
      call plts(sigmai, sigmaf, pi, pf, kc+1, title, 0.2, 1.2, 
     2          pmin, labx, pltloc)
c
      labx = 'pv/(2*omega)'
      pltloc = 'ur'
      write(title,230) lat
  230 format('pv/(2*omega) at ',f4.0)
      call plts(pvi, pvf, pi, pf, kc+1, title, 0., 1., 
     2          pmin, labx, pltloc)
c
      labx = 'zeta/(2*omega)'
      pltloc = 'lr'
      write(title,240) lat
  240 format('zetar/(2*omega) (x10**5) at ',f4.0)
      call plts(zetri, zetrf, pi, pf, kc+1, title,-10.,10., 
     2          pmin, labx, pltloc)
      go to 10
  999 continue
c
c   close workstation window
c
      call clsgks
c
      end
      subroutine plts(fi, ff, pi, pf, n, title, xmin, xmax, pmin,
     2                labx, pltloc )
      integer                         n
      real           fi(*), ff(*), pi(*), pf(*)                     
      character   labx*20, laby*20, title*(*), pltloc*2
      real pmaj(5)
      data pmaj/1000., 900., 800., 700., 600./ 
      print *, 'plotting: ',title
c
c   establish dash pattern for plots
c
      if(pltloc .eq. 'ur') then
         fb = 0.55
         ft = 0.90
         fl = 0.55
         fr = 0.95
      elseif(pltloc .eq. 'ul') then
         fb = 0.55
         ft = 0.90
         fl = 0.10
         fr = 0.50
      elseif(pltloc .eq. 'lr') then
         fb = 0.05
         ft = 0.40
         fl = 0.55
         fr = 0.95
      elseif(pltloc .eq. 'll') then
         fb = 0.05
         ft = 0.40
         fl = 0.10
         fr = 0.50
      endif
      pmn = pmin
      pmx = pi(1) 
      if(xmin.eq.0. .and. xmax.eq.0.) then
         call mima(n, ff, xmin, xmax, -99.99)
      endif
      print *, xmin, xmax, pmn, pmx   
      xl = xmin 
      xr = xmax
      yb = pmx
      yt = pmn  
      call set(fl, fr, fb, ft, xl, xr, yb, yt, 2)
      call pfseti('ISL', -3)
      call pfseti('ISR', -1)
      call pfseti('ISB',  2)
      laby = 'Pressure (mb)'
      call ptitle(title, 1.0)
      call perimn(labx, laby)
      if (pltloc .eq. 'ur' .or. pltloc .eq. 'lr') then
         call yscalu(1, 1, laby, 5, pmaj, pmaj, 0, smin, 0)
      else
         call yscalu(3, 1, laby, 5, pmaj, pmaj, 0, smin, 0)
      endif
      isold = ior (ishift(32767,1), 1)
      idash = ior( ishift(2184, 1), 1)
      do 33 i=1,n
         if(pi(i) .gt. pmin) then
            nn = i 
         endif
   33 continue
      call dashdb(isold)
      call curved(fi, pi, nn)
      call dashdb(idash)
      call curved(ff, pf, nn)
      if(pltloc .eq. 'lr') call frame
      return
      end
      SUBROUTINE MIMA( N, X, XMIN, XMAX, XMIS )
      INTEGER            N
      REAL               X(N), XMIN, XMAX, xmis
C
C   Finds the minimum and maximum values in the vector  X  of length  N
C
      INTEGER    I
      REAL       BIG
      PARAMETER  ( BIG = 3.0E38 )
C
      XMIN = +BIG
      XMAX = -BIG
      DO 10 I=1,N
          if(x(I) .eq. xmis) go to 10
          XMIN = MIN( XMIN, X(I) )
          XMAX = MAX( XMAX, X(I) )
   10 CONTINUE
      RETURN
      END
