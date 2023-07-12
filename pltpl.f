      subroutine pltpl(u, p, jc2, jtot, kc, pmin, pmax, 
     2                 tmin, tmax, time, jmin, jplt)     
c
c   plot out fields in theta, potential latitude space
c
      integer  jmin, jplt, jtot, i
      integer  jc2, kc, madx, mady, midx, midy, nsdx, nsdy
c
      real  u(-jc2:jc2,0:kc), p(-jc2:jc2,0:kc)
      real  pmin, pmax, tmin, tmax, convr, cl, finc, fmin
      character time*(*)
      integer iara(10), igra(10), iama(40000)
      integer nmaj, nmin, ltype
      integer iasf(13)
      parameter (nmaj = 5, nmin=8)
      real xmaj(nmaj), smaj(nmaj), xmin(nmin), smin(nmin)
      real rwrk(5000), xcra(1000), ycra(1000)
      integer iwrk(1000)
      real x1, x2, y1, y2, fr, fl, ft, fb
c
      character labx*80, laby*80, title*80
      external shader, drawcl
      data    xmaj/30., 15., 0., -15., -30./
      data    xmin/25., 20., 10., 5., -5., -10., -20., -25./
      data    iasf/13*1/
c
      convr = acos(-1.) / 180.
      do 10 i=1,nmaj
         smaj(i) = sin(xmaj(i)*convr)
   10 continue
      do 12 i=1,nmin
         smin(i) = sin(xmin(i)*convr)
   12 continue
      call gsasf(iasf)
c
c   define labels and title for plot
c
      labx = 'potential latitude'
      laby = 'theta (K)'
      mady = 4
      midy = 2
      madx = 6
      midx = 2
      nsdx = 0
      nsdy = 0
      if(time .eq. 'initial') then
         write(title,60) time
   60    format('p(kPa) at ',a7,' time')
      else
         write(title,61) time
   61    format('p(kPa) and u(m/s) at ',a7,' time')
      endif
      call pfhq
      call pfdupl
      call pfsetr( 'RATIO', 0.6)  
      call setbig(pmin, pmax, tmin, tmax, 1)
      call pfseti('ISB', -3)
      call pfseti('IST', -1)
      call perimf(labx,laby,madx,mady,midx,midy,nsdx,nsdy)
      call xscalu(3, 1, labx, nmaj, smaj,xmaj,nmin,smin,xmin)
      call ptitle(title, 0.8)
      call getset(fl, fr, fb, ft, x1, x2, y1, y2, ltype)
      call cpseti('SET', 0)
      call cpsetr('XC1', x1)
      call cpsetr('XCM', x2)
      call cpsetr('YC1', y1)
      call cpsetr('YCN', y2)
      call cpseti('LIS', 4)
      call cpseti('RWC', 500)
      call pcseti('QU  - quality flag for labels', 0)
      call cpsetr('CIS - contour interval', 5.)
      call cprect(p(jmin,0),jtot,jplt,kc+1,rwrk,5000,iwrk,1000)
      call cpcldr(p(jmin,0),rwrk,iwrk)
      if(time .eq. 'initial') go to 999
c
c  contour winds
c
      fmin = -5.0
      finc = 0.5
      call cpseti('CLS - contour level selection flag', 0) 
      call cpseti('NCL - number of contour levels', 21) 
      do 100  i=1,21
         call cpseti('PAI - param idex', i)
         cl = fmin + finc*float(i-1)
         call cpsetr('CLV - contour level value', cl)
         if(mod(cl,finc) .eq. 0.0) then
            call cpseti('CLU - contour level use', 3)
         else
            call cpseti('CLU - contour level use', 1)
         endif
         if(cl .ne. 0.0) then
            call cpseti('AIA - area id above line', 0)
            call cpseti('AIB - area id below line', 0)
         else
            call cpseti('AIA - area id above line', 2)
            call cpseti('AIB - area id below line', 1)
         endif
         if(cl .lt. 0.) then                
            call cpseti('CLD - dash pattern', 61680)
         endif
  100 continue
      call cprect(u(jmin,0),jtot,jplt,kc+1,rwrk,5000,iwrk,1000)
      call arinam(iama, 40000)
c
c  don't plot highs and lows
c
      call cpseti('LLP - line label position', 3)
      call cpsetc('HIT - high char. string', ' ')
      call cpsetc('LOT - low  char. string', ' ')
      call cpsetr('ILY - info label y-coor', 1.2)
      call cplbam(u(jmin,0),rwrk,iwrk,iama)
      call cpcldm(u(jmin,0),rwrk,iwrk,iama,drawcl)
      call cplbdr(u(jmin,0),rwrk,iwrk)
c
c   add zero contour
c
      call cpclam(u(jmin,0),rwrk,iwrk,iama)
      call arscam(iama,xcra,ycra,1000,iara,igra,10,shader)
      call frame
  999 continue
      return
      end
