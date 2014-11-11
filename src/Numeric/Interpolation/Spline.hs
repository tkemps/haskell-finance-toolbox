{-# LANGUAGE ScopedTypeVariables,FlexibleContexts,StandaloneDeriving,UndecidableInstances #-}
module Numeric.Interpolation.Spline where

import Control.Monad
import Control.Monad.ST
import Data.Array.IArray
import Data.Array.MArray
import Data.Array.ST
import Data.STRef
import Numeric.Util

data (IArray a e,Fractional e) => Spline a e = 
  Spline {
    spline'x :: a Int e,
    spline'y :: a Int e,
    spline'y'' :: a Int e
    }

deriving instance (IArray a e,Fractional e,Show (a Int e)) => Show (Spline a e)

spline :: forall a e.(Fractional e,IArray a (e,e),IArray a e)
         => a Int (e,e) -> Spline a e
spline v = runST $ (do
  y'' <- newArray (l,n) 0 :: ST s (STArray s Int e) -- second derivatives
  u <- newArray (l,n-1) 0 :: ST s (STArray s Int e) -- temporary storage
  forM_ [l+1..n-1]
    (\i->do
          let xip1 = x (i+1)
          let xi   = x i
          let xim1 = x (i-1)
          let yip1 = y (i+1)
          let yi   = y i
          let yim1 = y (i-1)
          let sig = (xi - xim1) / (xip1 - xim1)
          y''im1 <- readArray y'' (i-1)
          let p = sig *  y''im1 + 2.0
          writeArray y'' i ((sig - 1) / p)
          let z = (yip1 - yi) / (xip1 - xi) - (yi - yim1) / (xi - xim1)
          uim1 <- readArray u (i-1)
          writeArray u i ((6*z/(xip1 - xim1) - sig*uim1)/p))
  let qn = 0
  let un = 0
  y''nm1 <- readArray y'' (n-1)
  unm1 <- readArray u (n-1)
  writeArray y'' n ((un-qn*unm1)/(qn*y''nm1+1.0))
  forM_ (reverse [l..n-1])
    (\i->do
        y''i <- readArray y'' i
        y''ip1 <- readArray y'' (i+1)
        ui <- readArray u i
        let z = y''i*y''ip1+ui
        writeArray y'' i z)
  fy'' <- freeze y''
  return (Spline (amap fst v) (amap snd v) fy''))
  where (l,n) = bounds v
        x :: Int -> e
        x i = fst (v!i)
        y :: Int -> e
        y i = snd (v!i)

splineEval :: forall a e.(Ord e,Fractional e,IArray a e) => Spline a e -> e -> Maybe e
splineEval (Spline x y y'') x0 = 
  case binaryIntervalSearchIArray x x0 of
    Just i -> let h = x!i - x!(i+1)
                  a = (x!i - x0) / h
                  b = (x0 - x!(i+1)) / h
                  y0 = a * y!(i+1) + b*y!i
                      + ((a*a*a-a)*y''!(i+1)
                      +  (b*b*b-b)*y''!i)*h*h/6.0
             in Just y0
    Nothing -> Nothing
