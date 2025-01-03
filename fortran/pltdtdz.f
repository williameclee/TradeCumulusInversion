      program pltdtdz    
c
c   plots vertical profile of lapse rate as generated by invtci
c   code
c
      integer    kc    
      parameter (kc = 60)    
      real   ti(0:kc), pi(0:kc), thetai(0:kc), tki(0:kc), zi(0:kc)
      real   tf(0:kc), pf(0:kc), thetaf(0:kc), tkf(0:kc), zf(0:kc)
      real   dtdzi(0:kc), dtdzf(0:kc)
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
   10 continue
      read(2,*,end=999) lat
      do 30 k=0,kc
         read(2,*) kk, thetai(k), pi(k), ti(k), zi(k)   
   30 continue
      do 50  k=0,kc
         if(k .eq. 0) then
            dt = -3.*ti(k) + 4.*ti(k+1) - ti(k+2)
            dz = -3.*zi(k) + 4.*zi(k+1) - zi(k+2)
         elseif(k .eq. kc) then
            dt =  3.*ti(k) - 4.*ti(k-1) + ti(k-2)
            dz =  3.*zi(k) - 4.*zi(k-1) + zi(k-2)
         else
            dt = ti(k+1) - ti(k-1)
            dz = zi(k+1) - zi(k-1)
         endif
         dtdzi(k) = dt/dz * 1000.
   50 continue
      read(2,*) lat
      do 130 k=0,kc
         read(2,*) kk, thetaf(k), pf(k), tf(k), zf(k)   
  130 continue
      do 150  k=0,kc
         if(k .eq. 0) then
            dt = -3.*tf(k) + 4.*tf(k+1) - tf(k+2)
            dz = -3.*zf(k) + 4.*zf(k+1) - zf(k+2)
         elseif(k .eq. kc) then
            dt =  3.*tf(k) - 4.*tf(k-1) + tf(k-2)
            dz =  3.*zf(k) - 4.*zf(k-1) + zf(k-2)
         else
            dt = tf(k+1) - tf(k-1)
            dz = zf(k+1) - zf(k-1)
         endif
         dtdzf(k) = dt/dz * 1000. 
  150 continue
c
c   plot fields              
c
      labx = 'lapse rate (K/km)'
      write(title,210) lat
  210 format('results  at ',f4.0)
      call plts(dtdzi, dtdzf, pi, pf, kc+1, title, -10.,2.,   
     2          pmin, labx)
      go to 10
  999 continue
c
c   close workstation window
c
      call clsgks
c
      end
      subroutine plts(fi, ff, pi, pf, n, title, xmin, xmax, pmin,
     2                labx )
      integer                         n
      real           fi(*), ff(*), pi(*), pf(*)                     
      character   labx*20, laby*20, title*(*)
      real pmaj(9), smin(8)
      data pmaj/1000., 950., 900., 850., 800., 750., 700., 
     2                 650., 600./ 
      data smin/ 975., 925., 875., 825., 775., 725., 675., 625./

      print *, 'plotting: ',title
c
c   establish dash pattern for plots
c
         fb = 0.15
         ft = 0.90
         fl = 0.15
         fr = 0.50
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
      call pfdupl
      call pfhq
      call set(fl, fr, fb, ft, xl, xr, yb, yt, 2)
      call pfseti('ISL', -3)
      call pfseti('ISR', -1)
      laby = 'Pressure (mb)'
      call ptitle(title, 1.0)
c     call perime(labx, laby)
      call perimf(labx, laby, 6, 4, 2, 4, 0, 0)
      if (pltloc .eq. 'ur' .or. pltloc .eq. 'lr') then
         call yscalu(1, 1, laby, 9, pmaj, pmaj, 0, smin, 0)
      else
         call yscalu(3, 1, laby, 9, pmaj, pmaj, 8, smin, 0)
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
