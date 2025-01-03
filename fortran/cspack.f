      SUBROUTINE CSSET( N, X, F, IBC1, IBCN, IERR )
C
C
C   PURPOSE
C
C       SETS UP A CUBIC SPLINE TO INTERPOLATE SPECIFIED DATA
C
C   ARGUMENTS
C
C       DIMENSION   X(N), F(N,4)
C
C       INPUT
C
C           N       NUMBER OF DATA POINTS  (N.GE.2)
C
C           X       ABSCISSAS (IN INCREASING ORDER)
C
C           F       ON INPUT,  F(I,1)  CONTAINS THE ORDINATE
C                   CORRESPONDING TO  X(I)  (I=1, ..., N).
C                   WITH SOME BOUNDARY CONDITIONS, ADDITIONAL
C                   INFORMATION MUST BE SUPPLIED IN  F(1,2)
C                   AND/OR  F(N,2)  AS DESCRIBED BELOW.
C
C           IBC1,   INTEGERS SPECIFYING THE BOUNDARY CONDITIONS
C           IBCN    TO BE APPLIED AT  X(1)  AND  X(N), RESPECTIVELY.
C                   IBC1 (IBCN)  IS INTERPRETED AS FOLLOWS:
C                       VALUE    CONDITION
C                         1    FIRST DERIVATIVE AT ENDPOINT SUPPLIED
C                              IN  F(1,2)  (F(N,2))
C                         2    SECOND DERIVATIVE AT ENDPOINT SUPPLIED
C                              IN  F(1,2)  (F(N,2))
C                         3    'NOT-A-KNOT' CONDITION:  CONTINUITY
C                              OF THE THIRD DERIVATIVE IMPOSED AT
C                              X(2)  (X(N-1)) --REQUIRES  N.GE.3
C                         4    SLOPE AT ENDPOINT ESTIMATED BY FITTING
C                              A CUBIC POLYNOMIAL TO THE FIRST (LAST)
C                              FOUR DATA POINTS--REQUIRES  N.GE.4
C
C       OUTPUT
C
C           F       CONTAINS THE POLYNOMIAL COEFFICIENTS OF THE CUBIC
C                   SPLINE.  THE INPUT VALUES  F(I,1)  (I=1, ..., N)
C                   REMAIN UNCHANGED.
C
C           IERR    ERROR FLAG  (ZERO ON NORMAL RETURN)
C
C   ERROR CONDITIONS
C
C       IERR = 1    N  TOO SMALL FOR CHOICE OF BOUNDARY CONDITIONS
C       IERR = 2    ABSCISSAS NOT IN INCREASING ORDER
C
C   USE
C
C       FIRST CALL  CSSET  TO SET UP THE CUBIC SPLINE.  THEN USE
C       THE FUNCTION  CSVAL  TO EVALUATE THE SPLINE OR ANY OF ITS
C       DERIVATIVES AS DESIRED.
C
C   METHOD
C
C       A TRIDIAGONAL SYSTEM IS SET UP FROM THE REQUIREMENT OF
C       CONTINUITY OF THE SECOND DERIVATIVE AT THE INTERIOR
C       BREAKPOINTS PLUS THE BOUNDARY CONDITIONS.  THE SYSTEM IS
C       SOLVED BY GAUSSIAN ELIMINATION FOR THE FIRST DERIVATIVES
C       AT THE BREAKPOINTS, FROM WHICH THE COEFFICIENTS OF THE
C       CUBIC SPLINE AS A PIECEWISE POLYNOMIAL ARE OBTAINED.
C
C   HISTORY
C
C       WRITTEN BY SCOTT R. FULTON (FEBRUARY, 1982)
C       BASED ON THE ROUTINE  CUBSPL  BY CARL DEBOOR
C
C   REFERENCE
C
C       DEBOOR, CARL (1978):  A PRACTICAL GUIDE TO SPLINES.
C               SPRINGER-VERLAG, NEW YORK, 392 PP.
C
C
      INTEGER  N, IBC1, IBCN, IERR, I, JBC1, JBCN
      REAL  X(N), F(N,4), DD(3), FACT
C
C   ARGUMENT CHECKS
C
      IERR = 0
      JBC1 = IBC1
      JBCN = IBCN
      IF ( IBC1.LT.1 .OR. IBC1.GT.4 )  JBC1 = 3
      IF ( IBCN.LT.1 .OR. IBCN.GT.4 )  JBCN = 3
      IF ( N.LT.MAX0( 2 , JBC1 , JBCN ) )  THEN
          IERR = 1
          RETURN
      END IF
      DO 10 I=2,N
          IF ( X(I-1).GE.X(I) )  THEN
              IERR = 2
              RETURN
          END IF
   10 CONTINUE
C
C   COMPUTE  X  DIFFERENCES AND FIRST DIVIDED DIFFERENCES OF  F
C
      DO 20 I=2,N
          F(I,3) = X(I) - X(I-1)
          F(I,4) = (F(I,1) - F(I-1,1))/F(I,3)
   20 CONTINUE
C
C   SET UP FIRST EQUATION FROM THE BOUNDARY CONDITION AT  X(1)
C
      IF ( JBC1.EQ.1 )  THEN
          F(1,4) = 1.0
          F(1,3) = 0.0
      ELSE IF ( JBC1.EQ.2 )  THEN
          F(1,4) = 2.0
          F(1,3) = 1.0
          F(1,2) = 3.0*F(2,4) - F(2,3)*F(1,2)/2.0
      ELSE IF ( JBC1.EQ.3 )  THEN
          F(1,4) = F(3,3)
          F(1,3) = F(2,3) + F(3,3)
          F(1,2) = ( F(3,3)*(F(2,3) + 2.0*F(1,3))*F(2,4)
     2             + F(2,3)*F(2,3)*F(3,4) )/F(1,3)
      ELSE
          F(1,4) = 1.0
          F(1,3) = 0.0
          DD(2) = (F(3,4) - F(2,4))/(F(2,3) + F(3,3))
          DD(3) = (F(4,4) - F(3,4))/(F(3,3) + F(4,3))
          DD(3) = ( DD(3) - DD(2) )/(F(2,3) + F(3,3) + F(4,3))
          F(1,2) = F(2,4) - F(2,3)*(DD(2) - (F(2,3)+F(3,3))*DD(3))
      END IF
C
C   SET UP EQUATIONS FOR CONTINUITY OF SECOND DERIVATIVE AT INTERIOR
C   BREAKPOINTS  X  AND DO FORWARD PASS OF GAUSSIAN ELIMINATION
C
      DO 30 I=2,N-1
          FACT = -F(I+1,3)/F(I-1,4)
          F(I,2) = FACT*F(I-1,2) + 3.0*(F(I+1,3)*F(I,4)+F(I,3)*F(I+1,4))
          F(I,4) = FACT*F(I-1,3) + 2.0*(F(I,3) + F(I+1,3))
   30 CONTINUE
C
C   SET UP LAST EQUATION FROM THE BOUNDARY CONDITION AT  X(N)
C   AND COMPLETE FORWARD PASS OF GAUSSIAN ELIMINATION
C
      IF ( JBCN.EQ.2 )  THEN
          F(N,2) = 3.0*F(N,4) + F(N,3)*F(N,2)/2.0
          FACT   = -1.0/F(N-1,4)
          F(N,4) =  FACT*F(N-1,3) + 2.0
          F(N,2) = (FACT*F(N-1,2) + F(N,2))/F(N,4)
      ELSE IF ( JBCN.EQ.3 )  THEN
          FACT   =  F(N-1,3) + F(N,3)
          DD(1)  = (F(N-1,1) - F(N-2,1))/F(N-1,3)
          F(N,2) = ( F(N-1,3)*(F(N,3) + 2.0*FACT)*F(N,4)
     2             + F(N,3)*F(N,3)*DD(1) )/FACT
          FACT   = -FACT/F(N-1,4)
          F(N,4) = (FACT + 1.0)*F(N-1,3)
          F(N,2) = (FACT*F(N-1,2) + F(N,2))/F(N,4)
      ELSE IF ( JBCN.EQ.4 )  THEN
          DD(3) = (F(N-2,1) - F(N-3,1))/F(N-2,3)
          DD(2) = (F(N-1,1) - F(N-2,1))/F(N-1,3)
          DD(3) = ( DD(2) - DD(3))/(F(N-2,3) + F(N-1,3))
          DD(2) = (F(N,4) - DD(2))/(F(N-1,3) + F(N  ,3))
          DD(3) = ( DD(2) - DD(3))/(F(N-2,3) + F(N-1,3) + F(N,3))
          F(N,2) = F(N,4) + F(N,3)*(DD(2) + (F(N-1,3)+F(N,3))*DD(3))
      END IF
C
C   DO BACKWARD PASS OF GAUSSIAN ELIMINATION
C
      DO 40 I=N-1,1,-1
   40 F(I,2) = (F(I,2) - F(I,3)*F(I+1,2))/F(I,4)
C
C   COMPUTE POLYNOMIAL COEFFICIENTS OF THE CUBIC SPLINE
C
      DO 50 I=2,N
          FACT  =  F(I,3)
          DD(1) = (F(I,1) - F(I-1,1))/FACT
          DD(3) =  F(I-1,2) + F(I,2) - 2.0*DD(1)
          F(I-1,3) = 2.0*(DD(1) - F(I-1,2) - DD(3))/FACT
          F(I-1,4) = 6.0*DD(3)/(FACT*FACT)
   50 CONTINUE
      RETURN
      END
      FUNCTION CSVAL( XVAL, MD, N, X, F, ISW )
C
C
C   PURPOSE
C
C       EVALUATES THE CUBIC SPLINE SET UP BY THE SUBROUTINE  CSSET
C
C   ARGUMENTS
C
C       DIMENSION   X(N), F(N,4)
C
C       INPUT
C
C           XVAL    ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
C
C           MD      INTEGER SPECIFYING THE DERIVATIVE DESIRED:
C                       MD = 0    SPLINE VALUE
C                       MD = 1    FIRST  DERIVATIVE
C                       MD = 2    SECOND DERIVATIVE
C                       MD = 3    THIRD  DERIVATIVE
C
C           N,X,F   UNCHANGED FROM CALL TO  CSSET
C
C           ISW     INTEGER SWITCH TO HELP LOCATE THE INTERVAL
C                   CONTAINING  XVAL:
C                       ISW.GT.0    ASCENDING SEARCH
C                       ISW.EQ.0    BISECTION SEARCH
C                       ISW.LT.0    DESCENDING SEARCH
C                   IF  ISW = +1 (-1), THE SEARCH STARTS WITH THE
C                   FIRST (LAST) INTERVAL.  THUS WHEN INTERPOLATING
C                   AT AN INCREASING (DECREASING) SEQUENCE OF
C                   POINTS  XVAL = Y(J)  (J=1, ..., K), SET
C                   ISW = +J (-J)  FOR MAXIMUM EFFICIENCY.
C
C       OUTPUT
C
C           CSVAL   RETURNS THE REQUESTED SPLINE VALUE OR
C                   DERIVATIVE.  EXTRAPOLATION IS USED WHEN  XVAL
C                   IS OUTSIDE THE INTERVAL  (X(1),X(N)),  AND THE
C                   THIRD DERIVATIVE IS TAKEN TO BE CONTINUOUS
C                   FROM THE RIGHT AT THE BREAKPOINTS.
C
C   USE
C
C       FIRST CALL  CSSET  TO SET UP THE CUBIC SPLINE.  THEN USE
C       THE FUNCTION  CSVAL  TO EVALUATE THE SPLINE OR ANY OF ITS
C       DERIVATIVES AS DESIRED.
C
C   METHOD
C
C       AN INDEX  I  IS DETERMINED AS DESCRIBED ABOVE SUCH THAT
C       X(I).LE.XVAL.LT.X(I+1).  THE DESIRED VALUE IS THEN OBTAINED
C       USING THE POLYNOMIAL REPRESENTATION OF THE CUBIC SPLINE
C       ON THAT INTERVAL.
C
C   HISTORY
C
C       WRITTEN BY SCOTT R. FULTON (FEBRUARY, 1982)
C       BASED ON THE ROUTINE  PPVALU  BY CARL DEBOOR
C
C   REFERENCE
C
C       DEBOOR, CARL (1978):  A PRACTICAL GUIDE TO SPLINES.
C               SPRINGER-VERLAG, NEW YORK, 392 PP.
C
C
      INTEGER  MD, N, ISW, I, IH, IM, J, K
      REAL  CSVAL, XVAL, X(N), F(N,4), DX, FACT
      SAVE  I
      DATA  I / 1 /
C
C   RETURN ZERO FOR FOURTH- AND HIGHER-ORDER DERIVATIVES
C
      J = MAX0( MD+1 , 1 )
      IF ( J.GT.4 )  THEN
          CSVAL = 0.0
          RETURN
      END IF
C
C   FIND THE  X  INTERVAL CONTAINING  XVAL
C
      IF ( ISW )  10, 30, 80
C
C   1. DESCENDING SEARCH
C
   10 I = MAX( I, 1 )
      I = MIN( I, N-1 )
      IF ( ISW.EQ.-1 .OR. XVAL.GE.X(I+1) )  I = N-1
      I = I+1
   20 I = I-1
      IF ( X(I).GT.XVAL .AND. I.GT.1 )  GO TO 20
      GO TO 100
C
C   2. BISECTION SEARCH
C
   30 I  = 1
      IH = N
   40 IF ( IH-I.LE.1 )  GO TO 100
      IM = (I + IH)/2
      IF ( XVAL - X(IM) )  50, 60, 70
   50 IH = IM
      GO TO 40
   60 I  = IM
      IH = IM
      GO TO 40
   70 I = IM
      GO TO 40
C
C   3. ASCENDING SEARCH
C
   80 I = MAX( I, 1 )
      I = MIN( I, N-1 )
      IF ( ISW.EQ.+1 .OR. X(I).GT.XVAL )  I = 1
      I = I-1
   90 I = I+1
      IF ( XVAL.GE.X(I+1) .AND. I+1.LT.N )  GO TO 90
C
C   EVALUATE THE CUBIC SPLINE
C
  100 CSVAL = F(I,4)
      IF ( J.EQ.4 )  RETURN
      FACT = 4-J
      DX = XVAL - X(I)
      DO 200 K=3,J,-1
          CSVAL = F(I,K) + DX*CSVAL/FACT
          FACT  = FACT - 1.0
  200 CONTINUE
      RETURN
      END
 

