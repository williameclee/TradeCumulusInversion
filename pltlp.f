      subroutine pltlp(t, jl, ne, pmin, pmax, pb, pt, time)
      integer  jl, ne, madx, mady, midx, midy, nsdx, nsdy, jtot 
c
c   routine to theta in pressure and latitude space
c
      real  t(-jl:jl,0:ne)
      real  pmin, pmax, pb, pt    
      character labx*80, laby*80, title*80, time*(*)
c
c   define labels and title for plot
c
      labx = 'latitude'
      laby = 'p (kPa)' 
      madx = 6
      midx = 2
c     call pfsetr( 'RATIO', 0.4)  
      call pfsetr( 'RATIO', 0.8)  
      mady = 4
      midy = 2
      nsdx = 0
      nsdy = 0
      jtot = 2*jl + 1
      write(title,65) time
   65 format('potential temperature (K) at ',a7,' time')
      call pfhq
      call pfdupl
      call setbig(pmin, pmax, pb, pt, 1)
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
      call cpcnrc(t, jtot, jtot, ne+1, 0.,0.,0.50,-1,-1,-924)
      call pfsetr( 'RATIO', 0.4)  
      call frame
      return
      end
