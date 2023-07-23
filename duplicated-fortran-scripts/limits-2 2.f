      Subroutine limits(scmin, scmax, sc, jc, convr, jmin, jmax, jplt)
c
c   compute grid point limits (i.e., jmin, jmax) for potential latitude
c   plots
c
      double precision sc(-jc:jc), scmin, scmax, convr, snmin, dif
      integer   jmin, jmax, jplt, jc, j
      snmin = 2.
      do 40 j=-jc,jc
         dif = abs(scmin-sc(j))
         if(dif .lt. snmin) then
            jmin = j/2
            snmin = dif
         end if
   40 continue
      snmin = 2.
      do 50 j=-jc,jc
         dif = abs(scmax-sc(j))
         if(dif .lt. snmin) then
            jmax = j/2
            snmin = dif
         end if
   50 continue
      jplt = jmax - jmin + 1
      return
      end
