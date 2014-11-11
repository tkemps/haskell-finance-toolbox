{-# OPTIONS_GHC -O0 #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-
C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE MACHAR(IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,
     1                   MAXEXP,EPS,EPSNEG,XMIN,XMAX)
C-----------------------------------------------------------------------
C  This Fortran 77 subroutine is intended to determine the parameters
C   of the floating-point arithmetic system specified below.  The
C   determination of the first three uses an extension of an algorithm
C   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
C   but not all, of the improvements suggested by M. Gentleman and S.
C   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
C   program was published in the book Software Manual for the
C   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
C   Englewood Cliffs, NJ, 1980.
C
C  The program as given here must be modified before compiling.  If
C   a single (double) precision version is desired, change all
C   occurrences of CS (CD) in columns 1 and 2 to blanks.
C
C  Parameter values reported are as follows:
C
C       IBETA   - the radix for the floating-point representation
C       IT      - the number of base IBETA digits in the floating-point
C                 significand
C       IRND    - 0 if floating-point addition chops
C                 1 if floating-point addition rounds, but not in the
C                   IEEE style
C                 2 if floating-point addition rounds in the IEEE style
C                 3 if floating-point addition chops, and there is
C                   partial underflow
C                 4 if floating-point addition rounds, but not in the
C                   IEEE style, and there is partial underflow
C                 5 if floating-point addition rounds in the IEEE style,
C                   and there is partial underflow
C       NGRD    - the number of guard digits for multiplication with
C                 truncating arithmetic.  It is
C                 0 if floating-point arithmetic rounds, or if it
C                   truncates and only  IT  base  IBETA digits
C                   participate in the post-normalization shift of the
C                   floating-point significand in multiplication;
C                 1 if floating-point arithmetic truncates and more
C                   than  IT  base  IBETA  digits participate in the
C                   post-normalization shift of the floating-point
C                   significand in multiplication.
C       MACHEP  - the largest negative integer such that
C                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that
C                 MACHEP is bounded below by  -(IT+3)
C       NEGEPS  - the largest negative integer such that
C                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
C                 NEGEPS is bounded below by  -(IT+3)
C       IEXP    - the number of bits (decimal places if IBETA = 10)
C                 reserved for the representation of the exponent
C                 (including the bias or sign) of a floating-point
C                 number
C       MINEXP  - the largest in magnitude negative integer such that
C                 FLOAT(IBETA)**MINEXP is positive and normalized
C       MAXEXP  - the smallest positive power of  BETA  that overflows
C       EPS     - the smallest positive floating-point number such
C                 that  1.0+EPS .NE. 1.0. In particular, if either
C                 IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
C                 Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2
C       EPSNEG  - A small positive floating-point number such that
C                 1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
C                 or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEPS.
C                 Otherwise,  EPSNEG = (IBETA**NEGEPS)/2.  Because
C                 NEGEPS is bounded below by -(IT+3), EPSNEG may not
C                 be the smallest number that can alter 1.0 by
C                 subtraction.
C       XMIN    - the smallest non-vanishing normalized floating-point
C                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
C       XMAX    - the largest finite floating-point number.  In
C                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
C                 Note - on some machines  XMAX  will be only the
C                 second, or perhaps third, largest number, being
C                 too small by 1 or 2 units in the last digit of
C                 the significand.
C
C     Latest revision - April 20, 1987
C
C     Author - W. J. Cody
C              Argonne National Laboratory
C
C-----------------------------------------------------------------------
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
CS    REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
CS   1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
CD    DOUBLE PRECISION A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
CD   1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CD    CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
  520 RETURN
C---------- LAST CARD OF MACHAR ----------
      END
C
C   TEST DRIVER FOR MACHAR
C
      INTEGER IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP
CS    REAL EPS,EPSNEG,XMIN,XMAX
CD    DOUBLE PRECISION EPS,EPSNEG,XMIN,XMAX
C
      WRITE(6,1001)
      CALL MACHAR(IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,MAXEXP,
     1      EPS,EPSNEG,XMIN,XMAX)
      WRITE (6,1000) IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,
     1   MAXEXP,EPS,EPSNEG,XMIN,XMAX
      STOP
 1000 FORMAT(19H OUTPUT FROM MACHAR//7H BETA =,I5/7H    T =,I5/
     1   7H  RND =,I5/7H NGRD =,I5/9H MACHEP =,I5/8H NEGEP =,I5/
     2   7H IEXP =,I5/9H MINEXP =,I5/9H MAXEXP =,I5/6H EPS =,
     3   E25.13/9H EPSNEG =,E25.13/7H XMIN =,E25.13/
     4   7H XMAX =,E25.13)
 1001 FORMAT(18H1GOING INTO MACHAR/)
C  ---------- LAST CARD OF DRIVER ----------
      END
-}
module Numeric.Machine
        (beta,
        mantissaLength,
        RoundingMethod(IEEE,Truncate,OtherRounding),roundingMethod,
        negep,
        epsneg,
        machep,
        eps,
        guardDigits,
        iexp,
        xmin,xmax,
        minexp,maxexp,
        partialUnderflow)
    where

import Data.Bits

one :: Fractional e => e
one = 1

two :: Fractional e => e
two = one+one

zero :: Fractional e => e
zero = one-one

-- |Smallest integral number a with rounding errors in (((a+a)+1)-(a+a))-1
malcolmsA :: forall e.(Eq e,Fractional e) => e
malcolmsA = f (one :: e)
    where f a = let a' = a+a
                    temp = a+one
                    temp1 = temp-a
                in if temp1-one==zero then f a' else a

-- |The function beta returns the radix for the floating-point representation.
beta :: forall a.(Fractional a,RealFrac a) => a -> Int
beta _ = f one
    where m :: a
          m = malcolmsA :: a
          f :: a -> Int
          f b = let b' = b + b
                    t = m + b'
                    i = truncate (t - malcolmsA)
                in if i==0 then f b' else i

betaInv :: forall e.(Fractional e,RealFrac e) => e
betaInv = 1/(fromRational.toRational $ beta (0::e))

mantissaLength :: forall a.(Fractional a,RealFrac a) => a -> Int
mantissaLength _ = f 0 one
    where f it b = if temp1-one==zero then f (it+1) b' else it+1
                where b' = b*(fromRational.toRational $ beta (0::a))
                      temp = b'+one
                      temp1 = temp-b'

data RoundingMethod = IEEE | Truncate | OtherRounding deriving (Eq,Show)

roundingMethod :: forall a.(Fractional a,RealFrac a) => a -> RoundingMethod
roundingMethod _ = case irnd of
                     0 -> Truncate
                     1 -> OtherRounding
                     2 -> IEEE
    where bet = (fromRational.toRational $ beta (0::a))
          betah = bet/two
          a = malcolmsA
          temp = a+betah
          irnd = if temp-a /= zero then 1 
                 else let tempa = a+bet
                          temp = tempa+betah
                      in if temp-tempa/=zero then 2
                         else 0

negepEpsneg :: forall e.(Fractional e,RealFrac e) => (Int,e,e)
negepEpsneg = (-n',e,a)
    where n = (mantissaLength (0::e))+3
          betaInv = one/(fromRational.toRational $ beta (0::e))
          a = (iterate (*betaInv) one)!!n
          f a n = if temp-one == zero then f a' n' else (n,a)
                where temp = one-a
                      a' = a*(fromRational.toRational $ beta (0::e))
                      n' = n-1
          (n',e) = f a n

first (x,_,_) = x
second (_,x,_) = x
third (_,_,x) = x

negep :: forall e.(Fractional e,RealFrac e) => e -> Int
negep _ = let x = negepEpsneg :: (Int,e,e)
          in first x

epsneg :: (Fractional e,RealFrac e) => e
epsneg = second negepEpsneg

machepEps :: forall e.(Fractional e,RealFrac e) => (Int,e)
machepEps = f (third negepEpsneg) (-(mantissaLength (0::e))-3)
    where f a m = if temp-one == zero then f a' m' else (m,a)
            where temp = one+a
                  a' = a*(fromRational.toRational $ beta (0::e))
                  m' = m+1

machep :: forall e.(Fractional e,RealFrac e) => e -> Int
machep _ = let x = machepEps :: (Int,e)
           in fst x

eps :: (Fractional e,RealFrac e) => e
eps = snd machepEps

guardDigits :: forall e.(Fractional e,RealFrac e) => e -> Int
guardDigits _ = let temp :: e = one + eps
                in if roundingMethod (0::e)==Truncate && temp*one-one /= zero then 1 else 0

betaFract :: Fractional e => e
betaFract = fromRational.toRational $ beta 0

f1 :: forall a.(RealFrac a) => a -> (Int,Int,a)
f1 _ = f1' 0 1 betaInv
    where t = one+eps
          f1' i k z = let z' = z*z
                          a = z'*one
                          temp = z'*t*betaInv
                      in if a+a==zero || (abs z')>z || temp*betaFract==z' -- underflow
                         then (i,k,z')
                         else f1' (i+1) (k+k) z'

f2 :: forall a t.(RealFrac a) => a -> Int -> Int -> (Int,Int)
f2 _ i k = if bet/=10 then (i+1,k+k)   -- returns (iexp,mx)
           else (iexp,iz+iz-1)
    where bet = beta (0::a)
          g iz iexp = if k>=iz then g (iz*bet) (iexp+1) else (iz,iexp)
          (iz,iexp) = g bet 2

other :: forall e.(Fractional e,RealFrac e,Ord e) => (e,e,Int,Int,Int,Bool)
other = (xmin, xmax''', iexp, minexp, maxexp''', partialUnderflow)
    where bet = beta (0::e)
          betFract = fromRational.toRational $ bet
          (i,k,y) = f1 (0::e)
          (iexp,mx) = f2 (0::e) i k
          f3 y k = if a+a/=zero && (abs y')<y  -- returns (k,xmin,partialUnderflow,a,y)
                       then if y'*(one+eps)*betaInv*betFract==y' && y'*(one+eps)/=y'
                               then (k+1,y',True,a,y')
                               else f3 y' (k+1)
                       else (k,y,False,a,y')
                where y' = y*betaInv
                      a = y'*one
          (k',xmin,partialUnderflow,a,y') = f3 y k
          minexp = -k'
          (mx',iexp') = if mx<=k+k-3 && bet/=10 then (mx+mx,iexp+1) else (mx,iexp)
          maxexp = mx'+minexp-(if roundingMethod (0::e)==IEEE || partialUnderflow then 2 else 0)
          i' = maxexp+minexp
          maxexp' = if bet==2 && i'/=0 then maxexp-1 else maxexp
          maxexp'' = if i'>20 then maxexp'-1 else maxexp'
          maxexp''' = if a/=y' then maxexp''-2 else maxexp''
          xmax = one-epsneg         -- something like 0.9999999999999999
          xmax' = if xmax*one /= xmax then one-(fromRational.toRational $ bet)*epsneg else xmax
          xmax'' = xmax' / (xmin*betFract*betFract*betFract)
          i'' = maxexp'''+minexp+3
          xmax''' = if beta (0::e)==2
                    then (iterate (+xmax'') xmax'')!!i''
                    else (iterate (*betFract) xmax'')!!i''

xmin                = (\(x,_,_,_,_,_)->x) other
xmax                = (\(_,x,_,_,_,_)->x) other

-- |The number of bits (decimal places if IBETA = 10) reserved for the representation of the exponent
-- (including the bias or sign) of a floating-point number.
iexp                = (\(_,_,x,_,_,_)->x) other

minexp              = (\(_,_,_,x,_,_)->x) other
maxexp              = (\(_,_,_,_,x,_)->x) other
partialUnderflow    = (\(_,_,_,_,_,x)->x) other
