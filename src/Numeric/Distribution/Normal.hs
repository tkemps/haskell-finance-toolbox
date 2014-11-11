module Numeric.Distribution.Normal (
  pdfNormal,cdfNormal,quantileNormal,quantileNormal1
  ) where

-- |This function returns the value of the cumulative normal distribution. The implementation 
-- follows George Marsagila: Evaluating the Normal Distribution, Journal of Statistical Software
-- July 2004, Volume 11, Issue 4.
--
-- see @http://www.jstatsoft.org/v11/a05/paper@.
cdfNormal :: Double -> Double
cdfNormal x = phi' x 0 x 1
  where q = x*x
        phi' s t b i = if (s==t) 
                       then 0.5+s*exp(-0.5*q-0.91893853320467274178)
                       else let i'=i+2
                                b'=b*q/i'
                                s'=s+b'
                            in phi' s' s b' i'

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
