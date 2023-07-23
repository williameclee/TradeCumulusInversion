      function erf( x )
      SINGLE_OR_DOUBLE erf, erfc, x
c
c   returns the value of the error function  erf(x)
c
c   uses rational approximation 7.1.26 of Abramowitz and Stegun
c   claimed maximum error:  1.5e-07
c
c   written by Scott R. Fulton  (07/20/87)
c   modified to avoid underflow (06/18/90)
c
      SINGLE_OR_DOUBLE  p, a1, a2, a3, a4, a5, t
      save              p, a1, a2, a3, a4, a5
      data              p  /  0.3275911   /
      data              a1 /  0.254829592 /
      data              a2 / -0.284496736 /
      data              a3 /  1.421413741 /
      data              a4 / -1.453152027 /
      data              a5 /  1.061405429 /
c
      if ( abs( x ).le.5.0 ) then
          t = 1.0/(1.0 + p*abs( x ))
          erf = 1.0 - t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))*exp( -x*x )
      else
          erf = 1.0
      end if
      if ( x.lt.0.0 )  erf = -erf
      return
c
c   returns the value of the complementary error function  erfc(x)
c
      entry erfc( x )
      if ( abs( x ).le.5.0 ) then
          t = 1.0/(1.0 + p*abs( x ))
          erfc = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*exp( -x*x )
      else
          erfc = 0.0
      end if
      if ( x.lt.0.0 )  erfc = 2.0 - erfc
      return
      end
