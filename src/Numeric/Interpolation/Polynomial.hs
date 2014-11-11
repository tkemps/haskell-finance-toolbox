module Numeric.Interpolation.Polynomial where

import Data.Array.IArray
import Data.List

-- |Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
-- an error estimate dy. If P(x) is the polynomial of degree N >= 1 such that P(xai) = yai;
-- i=1,...,n, then the returned value is y = P(x) together with an error estimate.
-- The list of points shall have the following form:
--      ps :: Array Int Point
--      ps = array (1,3) [(1,(1,2)),(2,(3,4)),(3,(7,1))]
polynomialInterpolation :: Array Int (Double,Double) -> Double -> (Double,Double)
polynomialInterpolation a x = (\l -> (sum l,last l)) $ (snd (a!ic)):step n0 (ic-1)
    where   (n0,n1) = bounds a
            closest a x = minimumBy (\(_,(x1,_)) (_,(x2,_)) -> compare (abs (x-x1)) (abs (x-x2))) (assocs a)
            ic = fst $ closest a x
            c m i = let xi = fst $ a!i
                        xipm = fst $ a!(i+m)
                    in if xi==xipm then error "polynomialInterpolation: Points must not have equal x-coordinates."
                       else if m==n0 then (xi-x)*((snd $ a!(i+1)) - (snd $ a!i))/(xi-xipm)
                       else (xi-x)*((c (m-1) (i+1)) - (d (m-1) i))/(xi-xipm)
            d m i = let xi = fst $ a!i
                        xipm = fst $ a!(i+m)
                    in if xi==xipm then error "polynomialInterpolation: Points must not have equal x-coordinates."
                       else if m==n0 then (xipm-x)*((snd $ a!(i+1)) - (snd $ a!i))/(xi-xipm)
                       else (xipm-x)*((c (m-1) (i+1)) - (d (m-1) i))/(xi-xipm)
            step m ic | m==n1-1 = if ic==0 then [c m 1] else [d m ic]   -- last step: c or d
                      | 2*(ic-n0)<n1-m = (c m (ic+1)):step (m+1) ic     -- c step
                      | otherwise = (d m ic):step (m+1) (ic-1)          -- d step
