module Numeric.Optimization.Brent1Dim where

import Data.Maybe
import Numeric.Machine

-- * One Dimensional Optimization
-- One dimensional unrestricted optimization of general functions that need not be linear or
-- differentiable.

-- |The function call @sign a b@ returns @(abs a)@ if @b>=0@ and @-(abs a)@ if @b<0@.
sign :: Double -> Double -> Double
sign a b = if b>=0
           then if a>=0 then a else -a
           else if a>=0 then -a else a

-- |Example: minimum1DimBrent 3 3.5 5 sin Nothing
minimum1DimBrent :: Double -> Double -> Double -> (Double -> Double) -> Maybe Double -> (Double, Double)
minimum1DimBrent a b c f tol = iterateBrent a' b' f (b,f b) (b,f b) (b,f b) Nothing Nothing tol' 1
    where a' = if a < c then a else c
          b' = if a > c then a else c
          tol' = fromMaybe eps tol

zEps = eps*1.0e-3
goldenRatio = (1+sqrt 5)/2
gold = 1/goldenRatio**2

-- |make a step towards the local minimum between a and b.
-- d is the last step width. e is the step width before d. Initially both are zero.
iterateBrent ::
  Double
  -> Double
  -> (Double -> Double)
  -> (Double, Double)
  -> (Double, Double)
  -> (Double, Double)
  -> Maybe Double
  -> Maybe Double
  -> Double
  -> Int
  -> (Double, Double)
iterateBrent a b f (x,fx) (v,fv) (w,fw) d e tol count =
        if count>countMax
        then error ("minimum1DimBrent: No conversion after "++show count++" iterations.")
        else if abs(x-0.5*(a+b))<=2*tol'-0.5*(b-a)
             then (x,fx)
             else let (u,fu,d',e') = improve a b f (x,fx) (v,fv) (w,fw) d e tol'
                  in if fu <= fx
                     then if u >= x
                          then iterateBrent x b f (u,fu) (w,fw) (x,fx) d' e' tol' (count+1)
                          else iterateBrent a x f (u,fu) (w,fw) (x,fx) d' e' tol' (count+1)
                     else if fu <= fw || w == x
                          then if u<x
                               then iterateBrent u b f (x,fx) (w,fw) (u,fu) d' e' tol' (count+1)
                               else iterateBrent a u f (x,fx) (w,fw) (u,fu) d' e' tol' (count+1)
                          else if u<x
                               then iterateBrent u b f (x,fx) (u,fu) (w,fw) d' e' tol' (count+1)
                               else iterateBrent a u f (x,fx) (u,fu) (w,fw) d' e' tol' (count+1)
    where countMax = 100
          tol' = tol * (abs x) + zEps

improve ::
  Double
  -> Double
  -> (Double -> t)
  -> (Double, Double)
  -> (Double, Double)
  -> (Double, Double)
  -> Maybe Double
  -> Maybe Double
  -> Double
  -> (Double, t, Maybe Double, Maybe Double)
improve a b f (x,fx) (v,fv) (w,fw) _ Nothing tol' = (u,f u,d',e')
    where (d',e') = makeGoldenSectionStep a b x
          u = makeStep x d' tol'
improve a b f (x,fx) (v,fv) (w,fw) d (Just e) tol' = (u,f u,d',e')
    where (d',e') = if (abs e>tol')
                  then tryParabolicStep a b (x,fx) (v,fv) (w,fw) d e tol'
                  else makeGoldenSectionStep a b x
          u = makeStep x d' tol'

makeStep :: Double -> Maybe Double -> Double -> Double
makeStep x Nothing tol' = x+(abs tol')
makeStep x (Just d) tol' = if abs d>=tol' then x+d else x+(sign tol' d)

tryParabolicStep ::
  Double
  -> Double
  -> (Double, Double)
  -> (Double, Double)
  -> (Double, Double)
  -> Maybe Double
  -> Double
  -> Double
  -> (Maybe Double, Maybe Double)
tryParabolicStep a b (x,fx) (v,fv) (w,fw) d e tol' =
        let r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (if (q > r) then -1 else 1)*((x-v)*q-(x-w)*r)
            q' = 2*abs(q-r)
        in if parabolicFitAcceptable a b x p q e
           then if x+p/q-a < 2*tol' || b-x-p/q < 2*tol'
                then (Just $ sign tol' (0.5*(a+b)-x),d)
                else (Just $ p/q,d)
           else makeGoldenSectionStep a b x
    where parabolicFitAcceptable a b x p q e = abs(p)<0.5*abs(q*e) && p>q*(a-x) && p<q*(b-x)

makeGoldenSectionStep :: Double -> Double -> Double -> (Maybe Double, Maybe Double)
makeGoldenSectionStep a b x =
        let e = if x >= 0.5*(a+b) then a-x else b-x        -- Here we take the golden section step into the larger of the two segments.
        in (Just $ gold*e,Just $ e)
