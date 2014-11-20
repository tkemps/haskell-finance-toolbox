module Numeric.Distribution.Normal (
  pdfNormal,cdfNormal,cdfNormal',quantileNormal,quantileNormal1,cdfBivariateNormal
  ) where

import           Data.Array.Repa ((!))
import           Data.Array.Repa (DIM1,DIM2)
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import           Data.Array.Repa.Shape

-- |This function returns the value of the cumulative normal distribution. The implementation 
-- follows George Marsagila: Evaluating the Normal Distribution, Journal of Statistical Software
-- July 2004, Volume 11, Issue 4.
--
-- see @http://www.jstatsoft.org/v11/a05/paper@.
cdfNormal' :: Double -> Double
cdfNormal' x = phi' x 0 x 1
  where q = x*x
        phi' s t b i = if (s==t) 
                       then 0.5+s*exp(-0.5*q-0.91893853320467274178)
                       else let i'=i+2
                                b'=b*q/i'
                                s'=s+b'
                            in phi' s' s b' i'

-- |Algorithm based on Hart (1968)
cdfNormal :: Double -> Double
cdfNormal x =
  let y = abs x
      exponential = exp $ -y^2 / 2
      f z = if x>0 then 1-z else z
  in if y>37
     then 0
     else if y<7.07106781186547
          then let sumA = (((((3.52624965998911E-02 * y
                               + 0.700383064443688) * y
                              + 6.37396220353165) * y
                             + 33.912866078383) * y
                            + 112.079291497871) * y
                           + 221.213596169931) * y 
                          + 220.206867912376
                   sumB = ((((((8.83883476483184E-02 * y
                                + 1.75566716318264) * y
                               + 16.064177579207) * y
                              + 86.7807322029461) * y
                             + 296.564248779674) * y
                            + 637.333633378831) * y
                           + 793.826512519948) * y 
                          + 440.413735824752
               in f $ exponential * sumA / sumB
          else let sumA = y + 1 /(y + 2 /(y + 3 /(y + 4 /(y + 0.65))))
               in f $ exponential / (sumA * 2.506628274631)

evlpl :: [Double] -> Double -> Double
evlpl [] x = 0
evlpl (a0:[]) x = a0
evlpl (a0:as) x = a0+x*(evlpl as x)

-- | The call 'quantileNormal p' returns x  such that cdfNormal x = p. We use the Newton
-- algorithm in order to search for the inverse of cdfNormal. The starting point is given
-- by 'quantileNormal1 p'.
--
-- See @http://www3.sympatico.ca/craymer/software/fortran/dcdflib/@.
quantileNormal :: Double      -- ^ The probability whose normal deviate is sought.
                  -> Double
quantileNormal p = if p<0 || p>1 then 0/0
                   else if p==0 then -(1/0)
                   else if p==1 then 1/0
                   else let pp = max (1-p) p
                            x0 = quantileNormal1 pp
                            newton _ 0 = x0 -- error "quantileNormal did not converge."
                            newton x n = let cum = cdfNormal x 
                                             dx = (cum-pp)/(pdfNormal x)
                                             x' = x-dx
                                         in if abs (dx/x')<eps then x' else newton x' (n-1)
                            x = newton x0 maxIter
                        in if p>=0.5 then x else -x
  where maxIter = 1000
        eps = 1e-13

-- | Quick approximation to normal quantile according to rational function on page 95 
-- of Kennedy and Gentle, Statistical Computing, Marcel Dekker, NY, 1980.
quantileNormal1 :: Double -> Double
quantileNormal1 p = if p<0 || p>1 then 0/0
                   else if p==0 then -(1/0)
                   else if p==1 then 1/0
                   else let (sign,z) = if p>0.5 then (-1,p) else (1,1-p)
                            y = sqrt(-2*log z)
                        in sign*(y + f y)
  where f x = let x2 = x*x
                  x3 = x2*x
                  x4 = x3*x
                  n = -0.322232431088 - x -0.342242088547*x2 -0.204231210245e-1*x3 
                      -0.453642210148e-4*x4
                  d = 0.993484626060e-1 + 0.588581570495*x + 0.531103462366*x2
                       + 0.103537752850*x3 + 0.38560700634e-2*x4
              in n/d

pdfNormal x = r2pi*exp (nhalf*x*x)
  where r2pi = 0.3989422804014326
        nhalf = -0.5

-- |This function is based on the method described by Drezner, Z and G.O. Wesolowsky, (1990),
-- On the computation of the bivariate normal integral, Journal of Statist. Comput. Simul.
-- 35, pp. 101-107, with major modifications for double precision, and for |rho| close to 1.

cdfBivariateNormal x y rho =
  let w = R.transpose $ R.fromListUnboxed (ix2 3 10) $ [
        0.17132449237917,  -- 0x0
        0.360761573048138, -- 1x0
        0.46791393457269,  -- 2x0
        0,
        0,
        0,
        0,
        0,
        0,
        0,                 -- 9x0
        
        4.71753363865118E-02,   -- 0x1
        0.106939325995318,      -- 1x1
        0.160078328543346,      -- 2x1
        0.203167426723066,
        0.233492536538355,
        0.249147045813403,
        0,
        0,
        0,
        0,                      -- 9x1
        
        1.76140071391521E-02,   -- 0x2
        4.06014298003869E-02,   -- 1x2
        6.26720483341091E-02,
        8.32767415767048E-02,
        0.10193011981724,
        0.118194531961518,
        0.131688638449177,
        0.142096109318382,
        0.149172986472604,
        0.152753387130726]      -- 9x2
      xx = R.transpose $ R.fromListUnboxed (ix2 3 10) $ [
        -0.932469514203152,
        -0.661209386466265,
        -0.238619186083197,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        
        -0.981560634246719,
        -0.904117256370475,
        -0.769902674194305,
        -0.587317954286617,
        -0.36783149899818,
        -0.125233408511469,
        0,
        0,
        0,
        0,
        
        -0.993128599185095,
        -0.963971927277914,
        -0.912234428251326,
        -0.839116971822219,
        -0.746331906460151,
        -0.636053680726515,
        -0.510867001950827,
        -0.37370608871542,
        -0.227785851141645,
        -7.65265211334973E-02]

      (ng,lg) = if abs rho<0.3
                then (1,3)
                else if abs rho<0.75
                     then (2,6)
                     else (3,10)

      h = -x
      k = -y
      hk = h * k
  in if abs rho<0.925
     then cdfNormal (-h)*cdfNormal (-k)
          + if abs rho>0
            then let hs = (h^2+k^2)/2
                     asr = asin rho
                     s = sum [let sn = sin (0.5*asr*(iss*xx!(ix2 i (ng-1))+1))
                              in w!(ix2 i (ng-1))*(exp ((sn*hk-hs)/(1-sn^2)))
                             | i<-[0..lg-1],iss<-[-1,1]]
                 in s*asr/(4*pi)
            else 0
     else let (k',hk') = if rho<0 then (-k,-hk) else (k,hk)
              bvn' = if abs rho<1
                     then let ass = (1-rho)*(1+rho)
                              a = sqrt ass
                              b = abs (h-k')
                              bs = b^2
                              c = (4-hk')/8
                              d = (12-hk')/16
                              asr = -0.5*(bs/ass+hk')
                              bvn0 = if asr > -100
                                     then a*(exp asr)
                                          *(1-c*(bs-ass)*(1-d*bs/5)/3
                                            + c*d*ass^2/5)
                                     else 0
                              bvn1 = bvn0 + if hk' > -100
                                            then -(exp (-hk'/2))
                                                 *(sqrt (2*pi))
                                                 *(cdfNormal (-b/a))
                                                 *b*(1-c*bs*(1-0.2*d*bs)/3)
                                            else 0
                              a' = a/2
                              s = sum [let xs = (a'*(iss*xx!(ix2 i (ng-1))+1))^2
                                           rs = sqrt (1-xs)
                                           asr = -(bs/xs+hk')/2
                                           ss = if asr > -100
                                                then a'*w!(ix2 i (ng-1))*(exp asr)
                                                     *((exp (-hk'*(1-rs)/(2*(1+rs))))/rs
                                                       - (1+c*xs*(1+d*xs)))
                                                else 0
                                       in ss
                                      | i<-[0..lg-1],iss<-[-1,1]]
                          in -(bvn1+s)/(2*pi)
                     else 0
          in if rho>0
              then bvn' + cdfNormal (-max h k')
              else -bvn' + if h<k'
                           then cdfNormal k' - cdfNormal h
                           else 0
{-
unsafePrint0 lbl x = unsafePerformIO $ do
  putStrLn (show x++"="++lbl)
  return 0.0
-}
