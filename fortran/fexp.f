      function fexp(x)
      real     fexp, x, k, fac
      data k /2.5608517/
c
      fac = exp(1./(x-1.))
      fexp = exp(-k*fac/x)
      return
      end
