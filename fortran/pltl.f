      subroutine pltlup(u, p, jl, kc, pmin, pmax, tmin, tmax, time)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy, jtot 
c
c   routine to plot zonal winds and pressure field in actual latitude, theta
c   space
c
      real  u(-jl:jl,0:kc), p(-jl:jl,0:kc)
      real  pmin, pmax, tmin, tmax
      character time*(*)
      character labx*80, laby*80, title*80
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'theta (K)'
      mady = 4
      midy = 2
      madx = 6
      midx = 1
      call pfsetr( 'RATIO', 0.4)  
c     madx = 4
c     midx = 3
c     call pfsetr( 'RATIO', 0.6)  
      nsdx = 0
      nsdy = 0
      jtot = 2*jl + 1
      call pfhq
      call pfdupl
      call setbig(pmin, pmax, tmin, tmax, 1)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call pcseti('QU  - quality flag for labels', 0)
      call cpseti('LLP - line label positioning', 3)
      call cpsetr('PC1 - penalty scheme constant 1', 1.5)
      call cpsetr('PC6 - penalty scheme constant 6', 1.2)
      call cpseti('LIS - line interval specifier ', 2)
c     call cpseti('LLS - line label size ', .016)
      call cpsetr('ILY - info label y-coor', 1.3)
      call cpseti('RWC - contour workspace', 500)
      call cpseti('T2D - tension on 2d smoother', 4)
      if(time .eq. 'initial') then 
         write(title,60) time
   60    format('p(kPa) ',a7,' time')
         call ptitle(title, 0.8)
         go to 70
      else
         write(title,61) time
   61    format('p(kPa) and u(m/s) at ',a7,' time')
         call ptitle(title, 0.8)
      endif
      call cpcnrc(u, jtot, jtot, kc+1, 0.,0., 0.25,-1,-1,-924)
c
c   plot pressure
c
   70 continue
      call cpsetr('ILY - info label y-coor', 1.4)
      call cpcnrc(p, jtot, jtot, kc+1, 0.,0., 5.,-1,-1,-924)
      call frame
      return
      end
      subroutine pltls(s, jl, kc, pmin, pmax, tmin, tmax, time)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy, jtot 
      integer  j, k
c
c   routine to plot pseudo-density
c
      real  s(-jl:jl,0:kc)
      real  pmin, pmax, tmin, tmax
      character labx*80, laby*80, title*80, time*(*)
      do 10 k=0,kc
      do 10 j=-jl,jl
         s(j,k)=s(j,k) + .01
   10 continue
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'theta (K)'
      madx = 6
      midx = 1
      call pfsetr( 'RATIO', 0.4)  
      mady = 3
      midy = 2
c     madx = 4
c     midx = 3
c     call pfsetr( 'RATIO', 0.6)  
      nsdx = 0
      nsdy = 0
      jtot = 2*jl + 1
      write(title,65) time
   65 format('normalized potential pseudodensity at ',
     2        a7,' time')
      call pfhq
      call pfdupl
      call setbig(pmin, pmax, tmin, tmax, 1)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call ptitle(title, 0.8)
      call pcseti('QU  - quality flag for labels', 0)
      call cpseti('LLP - line label positioning', 3)
      call cpsetr('PC1 - penalty scheme constant 1', 3.5)
      call cpsetr('PC3 - penalty scheme constant 3', 180.)
c     call cpsetr('PC4 - penalty scheme constant 4', .5)
      call cpsetr('PC6 - penalty scheme constant 6', 1.9)
      call cpseti('LIS - line interval specifier ', 2)
c     call cpseti('LLS - line label size ', .016)
      call cpsetr('ILY - info label y-coor', 1.3)
      call cpseti('RWC - contour workspace', 500)
      call cpseti('T2D - tension on 2d smoother', 4)
      call cpcnrc(s, jtot, jtot, kc+1, 0.,0., 0.2,-1,-1,-924)
      call frame
      return
      end
      subroutine pltlpv(pv, dpv, jl, kc, pmin, pmax, tmin, tmax, time)
      integer  ncl
      parameter (ncl = 5)
      integer  jl, kc, madx, mady, midx, midy, nsdx, nsdy, jtot 
c
c   routine to plot pv and shade area where its derivative is 
c   negative
c
      real  pv(-jl:jl,0:kc), dpv(-jl:jl,0:kc)
      real  pmin, pmax, tmin, tmax
      real rwrk(5000), xcra(1000), ycra(1000)
      real  fl, fr, fb, ft, x1, x2, y1, y2, cl
      integer iara(10), igra(10), iama(50000)
      integer ltype, i
      integer iasf(13)
      integer iwrk(1000)
      external shader, drawcl
      data    iasf/13*1/
      character labx*80, laby*80, title*80, time*(*)
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'theta (K)'
      mady = 3
      midy = 2
c     madx = 4
c     midx = 3
c     call pfsetr( 'RATIO', 0.6)  
      madx = 6
      midx = 1
      call pfsetr( 'RATIO', 0.4)  
      nsdx = 0
      nsdy = 0
      jtot = 2*jl + 1
      write(title,60) time
   60 format('normalized potential vorticity at ',a7,' time')
      call pfhq
      call pfdupl
      call setbig(pmin, pmax, tmin, tmax, 1)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call ptitle(title, 0.8)
      call pcseti('QU  - quality flag for labels', 0)
      call cpseti('LLP - line label positioning', 3)
      call cpsetr('PC1 - penalty scheme constant 1', 1.5)
c     call cpsetr('PC6 - penalty scheme constant 6', 1.2)
      call cpseti('LIS - line interval specifier ', 2)
c     call cpseti('LLS - line label size ', .016)
      call cpsetr('ILY - info label y-coor', 1.3)
      call cpseti('T2D - tension on 2d smoother', 4)
      call cpseti('RWC - contour workspace', 500)
c
c   plot pv
c
      call cpcnrc(pv, jtot, jtot, kc+1, 0.,0.,.1,-1,-1,-924)
c
c   shade areas of negative pv gradient 
c
      call gsasf(iasf)
      call getset(fl, fr, fb, ft, x1, x2, y1, y2, ltype)
      call cpseti('SET', 0)
      call cpsetr('XC1', x1)
      call cpsetr('XCM', x2)
      call cpsetr('YC1', y1)
      call cpsetr('YCN', y2)
      call cpseti('LIS', 4)
      call pcseti('QU  - quality flag for labels', 0)
      call cpsetr('CIS - contour interval', 5.)
      call cpseti('CLS - contour level selection flag', 0) 
      call cpseti('NCL - number of contour levels', ncl) 
      do 100  i=1,ncl
         call cpseti('PAI - param idex', i)
         cl = -1. + 0.5*float(i-1)
         call cpsetr('CLV - contour level value', cl)
         if(cl .ne. 0.0) then
            call cpseti('AIA - area id above line', 0)
            call cpseti('AIB - area id below line', 0)
         else
            call cpseti('AIA - area id above line', 2)
            call cpseti('AIB - area id below line', 1)
         endif
  100 continue
      call cprect(dpv,jtot,jtot,kc+1,rwrk,5000,iwrk,1000)
      call arinam(iama, 40000)
c
c  don't plot highs and lows
c
      call cpseti('LLP - line label position', 3)
      call cpsetc('HIT - high char. string', '')
      call cpsetc('LOT - low  char. string', '')
      call cpsetr('ILY - info label y-coor', 1.3)
c     call cplbam(dpv,rwrk,iwrk,iama)
c     call cpcldm(dpv,rwrk,iwrk,iama,drawcl)
c     call cplbdr(dpv,rwrk,iwrk)
c
c   add zero contour
c
      print *, 'calling routines to do shading'
      call cpclam(dpv,rwrk,iwrk,iama)
      call arscam(iama,xcra,ycra,1000,iara,igra,10,shader)
      call frame
c
c   contour plot of dpv
c
      write(title,65) time
   65 format('potential vorticity gradient at ',a7,' time')
      call pfhq
      call pfdupl
      call setbig(pmin, pmax, tmin, tmax, 1)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call ptitle(title, 0.8)
      call pcseti('QU  - quality flag for labels', 0)
      call cpseti('LLP - line label positioning', 3)
      call cpsetr('PC1 - penalty scheme constant 1', 1.5)
      call cpseti('LIS - line interval specifier ', 2)
c     call cpseti('LLS - line label size ', .016)
      call cpsetr('ILY - info label y-coor', 1.3)
      call cpseti('T2D - tension on 2d smoother', 4)
      call cpseti('RWC - contour workspace', 500)
c
c   plot dpv
c
      call cpcnrc(dpv, jtot, jtot, kc+1, 0.,0.,0.,-1,-1,-924)
      call frame
      return
      end
