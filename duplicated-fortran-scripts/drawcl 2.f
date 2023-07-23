      subroutine drawcl(xcs, ycs, ncs, iai, iag, nai)
      integer i, nai, iag(*), iai(*), idr
      integer ncs
      real    xcs(*), ycs(*)
      idr = 1
      do 101 i=1,nai
         if (iai(i).lt. 0) idr = 0
  101 continue
      if(idr.ne.0) call curved(xcs, ycs, ncs)
      return
      end
