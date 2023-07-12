      subroutine intdde( method, n, h, u0, f, u )
c
      integer  method, n
      real h, u0, f(0:n), u(0:n)
c
c=======================================================================
c intdde                                             Rick Taft, 4-1-1992
c
c This routine integrates the discrete ODE IVP
c
c          f(i) = du/dx  over some interval [a,b]
c          u(a) = u0
c
c to obtain the discrete function  u(i)  at the same discrete points as
c f(i)  is known.
c
c Inputs:
c ------
c     method = method of integration to use:
c                  1 = use cubic interpolation
c                  2 = use trapezoidal rule
c     n      = number of grid intervals between limits of integration
c                  (assumption:  n > 2)
c     h      = grid spacing = (b - a)/n
c     f(i)   = f(a + i*h)  for  i = 0,1,...,n
c     u0     = u(a)
c
c Output:
c ------
c     u(i) = u(a + i*h)  for  i = 0,1,...,n  such that  f = du/dx
c
c Note:
c ----
c     The problem is translated to the origin so the values  a,b  are
c     not needed.  The new domain is the interval [0,n*h].
c=======================================================================
c
c local variables
c
      integer  i
      real     hby2, hby24
c
c select method of integration
c
c========================== Cubic Interpolation ========================
c
      if ( method.eq.1 ) then
c
c perform first step to get u(1):
c   - fit cubic polynomial p(x) to f(0),f(1),f(2),f(3)
c   - then  u(1) = u(0) + integral of p(x) from 0 to h
c
          hby24 = h/24.0
          u(0) = u0
          u(1) = u0 + hby24*( 9.0*f(0) + 19.0*f(1) - 5.0*f(2) + f(3) )
c
c perform next step repeatedly to get u(i) for  i = 2,..,n-1
c   - fit cubic polynomial p(x) to f(i-2),f(i-1),f(i),f(i+1)
c   - then  u(i) = u(0) + integral of f from 0 to i*h
c                = u(i-1) + integral of p(x) from (i-1)*h to i*h
c
          do 100 i=2,n-1
	      u(i) = u(i-1)
     +             + hby24*( -f(i-2) + 13.0*( f(i-1) + f(i) ) - f(i+1) )
  100     continue
c
c perform last step to get u(n):
c   - fit cubic polynomial p(x) to f(n-3),f(n-2),f(n-1),f(n)
c   - then  u(n) = u(0) + integral of f from 0 to n*h
c                = u(n-1) + integral of p(x) from (n-1)*h to n*h
c
          u(n) = u(n-1)
     +         + hby24*( f(n-3) - 5.0*f(n-2) + 19.0*f(n-1) + 9.0*f(n) )
c
c========================== Trapezoidal Rule ===========================
c
      else
          hby2 = 0.5*h
          u(0) = u0
          do 200 i=1,n
              u(i) = u(i-1) + hby2*(f(i-1) + f(i))
  200     continue
      endif
c
      return
      end
