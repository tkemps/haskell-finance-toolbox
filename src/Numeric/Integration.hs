module Numeric.Integration (trapezoidal,qromb,qromb') where

import Control.Exception (throw,ArithException(LossOfPrecision))
import Control.Monad (forM_,when)
import Control.Monad.Cont (ContT,runContT,callCC,lift)
import Control.Monad.ST (ST,runST)
import Data.Array (Array)
import Data.Array.IArray (IArray,listArray,Ix)
import Data.Array.ST (STArray,MArray,newArray,newListArray,writeArray,
                      readArray,getElems,getBounds)
import Numeric.Interpolation.Polynomial (polynomialInterpolation)

trapezoidal :: (Double -> Double) -> Double -> Double -> [Double]
trapezoidal f a b = s1:(trapzd 2 s1)
  where s1 = 0.5*(b-a)*(f a + f b)
        trapzd n s = let d = (b-a)/(fromIntegral n)
                         sm = sum [f ((fromIntegral j + 0.5)*d+a) | j<-[1..n]]
                         s' = 0.5*(s+(b-a)*sm/(fromIntegral n))
                     in s':(trapzd (2*n) s')

-- |Returns the integral of the function f from a to b. Integration is performed 
-- by Romberg’s method of order 2k, where, e.g., k=2 is Simpson’s rule. Here 
-- k = 5 is used.
qromb :: (Double -> Double) -> Double -> Double -> Double
qromb f a b = runST $ runContT (qromb' f a b) return

qromb' :: (Double -> Double) -> Double -> Double -> ContT r (ST s) Double
qromb' f a b = callCC $ \exit -> do
        s <- lift (newArray (0,jMax) 0 :: ST s (STArray s Int Double))
        h <- lift (newArray (0,jMax+1) 0 :: ST s (STArray s Int Double))
        lift $ writeArray h 1 1.0
        forM_ [1..jMax] $ \j->do
          lift $ writeArray s j (trapzd!!j)
          when (j>=k) $ do
            ps <- lift $ zipArrayWithML (,) h s
            let v = listArray (0,k) (take k (drop (j-k) ps))
            -- extrapolate to limit h->0:
            let (ss,dss) = polynomialInterpolation v 0
            when (abs dss <= eps*abs ss) $ exit ss
          hjm1 <- lift $ readArray h (j-1)
          -- The following is a key step: The factor is 0.25 even  
          -- though the stepsize is decreased by only 0.5. This  
          -- makes the extrapolation a polynomial in h^2, not just 
          -- a polynomial in h.
          lift $ writeArray h j (0.25*hjm1)
        throw LossOfPrecision
  where eps = 1.0e-6   -- fractional accuracy desired
        jMax = 20      -- limits the total number of steps
        k = 5          -- number of points used in the extrapolation
        trapzd = trapezoidal f a b

{-# INLINE zipArrayWithMM #-}
-- | Constructs a new mutable array derived from the original arrays by applying a
-- function to each pair of the elements.
zipArrayWithMM :: (Ix i,MArray a1 e1 m,MArray a2 e2 m,MArray a3 e3 m) =>
                 (e1 -> e2 -> e3) -> a1 i e1 -> a2 i e2 -> m (a3 i e3)
zipArrayWithMM f a1 a2 = do 
  l1 <- getElems a1
  l2 <- getElems a2
  (n1,n2) <- getBounds a1
  newListArray (n1,n2) (zipWith f l1 l2)

{-# INLINE zipArrayWithMI #-}
-- | Constructs a new immutable array derived from the original arrays by applying a
-- function to each pair of the elements.
zipArrayWithMI :: (Ix i,MArray a1 e1 m,MArray a2 e2 m,IArray a3 e3) =>
                 (e1 -> e2 -> e3) -> a1 i e1 -> a2 i e2 -> m (a3 i e3)
zipArrayWithMI f a1 a2 = do
  (n1,n2) <- getBounds a1
  l <- zipArrayWithML f a1 a2
  return (listArray (n1,n2) l)

{-# INLINE zipArrayWithML #-}
-- | Constructs a list derived from the original arrays by applying a
-- function to each pair of the elements.
zipArrayWithML :: (Ix i,MArray a1 e1 m,MArray a2 e2 m) =>
                 (e1 -> e2 -> e3) -> a1 i e1 -> a2 i e2 -> m [e3]
zipArrayWithML f a1 a2 = do
  l1 <- getElems a1
  l2 <- getElems a2
  return (zipWith f l1 l2)
