      subroutine shader(xcs, ycs, ncs, iai, iag, nai)
      integer i, nai, iag(*), iai(*), ind(1200), ish
      integer ncs
      real    xcs(*), ycs(*), dst(1100)
      ish = 0
      do 101 i=1,nai
         if (iag(i).eq.3 .and. iai(i).eq.1) ish=1
  101 continue
      if(ish.ne.0) then
         call sfseti('ANGLE of fill lines', 45)
         call sfsetr('SPACING of fill lines', .004)
         call sfseti('DOT-FILL flag', 1)
         call sfwrld(xcs,ycs,ncs-1,dst,1100,ind,1200)
      endif
      return
      end
