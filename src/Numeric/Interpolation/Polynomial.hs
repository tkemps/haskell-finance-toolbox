{-# LANGUAGE FlexibleContexts #-}
module Numeric.Interpolation.Polynomial (
  linearInterpolation, logLinearInterpolation,cubicInterpolation,polynomialInterpolation
  ) where

import           Control.Applicative
import           Control.Monad.Identity
import           Control.Monad.ST
import qualified Data.Array.IArray as A
import           Data.Array.Repa ((!),DIM1,DIM2)
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import qualified Data.Array.Repa.Operators.Mapping as R
import           Data.Array.Repa.Shape
import           Data.Array.Repa.Slice
import           Data.List  (sortBy)
import qualified Data.List as L
import qualified Data.Vector.Unboxed as V hiding (length)
import qualified Data.Vector.Unboxed.Mutable as V
import           System.IO.Unsafe

type Vec s = R.Array s DIM1 Double
type VecD = Vec R.D
type VecU = Vec R.U

(!.) :: R.Source r e => R.Array r DIM1 e -> Int -> e
v !. i = v!ix1 i

(!..) :: R.Source r e => R.Array r DIM2 e -> (Int,Int) -> e
m !.. (i,j) = m!ix2 i j

vecFromList l = let n = length l
                    sh = R.shapeOfList [n] :: DIM1
                in R.fromListUnboxed sh l

brack :: (R.Source s Double) => Vec s -> Double -> (Int,Int)
brack x x0 = if n<=1
                then error "I need at least two points."
                else let ic = closest x x0
                     in if ic==0 then (0,1)
                        else if ic==n-1
                             then (n-2,n-1)
                             else if x0<=x!.ic
                                  then (ic-1,ic)
                                  else (ic,ic+1)
  where n=R.size (R.extent x)

brack4 :: (R.Source s Double) => Vec s -> Double -> (Int,Int)
brack4 x x0 = if n<=3
              then error "I need at least four points."
              else let (i1,i2) = brack x x0
                   in if i1==0
                      then (0,3)
                      else if i2==n-1
                           then (n-4,n-1)
                           else (i1-1,i2+1)
  where n=R.size (R.extent x)

cubicInterpolation :: (R.Source s Double) =>
  Vec s      -- ^Array of length at least 2 with x coordinates of known points
  -> Vec s    -- ^Array of length at least 2 with x coordinates of known points
  -> Double          -- ^x coordinate of new point
  -> Double
cubicInterpolation x y x0 = uncurry interp (brack4 x x0)
  where interp :: Int -> Int -> Double
        interp i0 i1 = sum [y!.i*product [(x0-x!.j)/(x!.i-x!.j) | j<-[i0..i1],i/=j] | i<-[i0..i1]]

linearInterpolation :: (R.Source s Double) =>
  Vec s      -- ^Array of length at least 2 with x coordinates of known points
  -> Vec s    -- ^Array of length at least 2 with x coordinates of known points
  -> Double          -- ^x coordinate of new point
  -> Double
linearInterpolation x y x0 = uncurry interp (brack x x0)
  where interp i0 i1 = (y!.i1-y!.i0) * (x0-x!.i0) / (x!.i1-x!.i0) + y!.i0

logLinearInterpolation :: (R.Source s Double) =>
  Vec s      -- ^Array of length at least 2 with x coordinates of known points
  -> Vec s    -- ^Array of length at least 2 with x coordinates of known points
  -> Double          -- ^x coordinate of new point
  -> Double
logLinearInterpolation x y x0 = uncurry interp (brack x x0)
  where interp i0 i1 = (y!.i1/y!.i0)**((x0-x!.i0)/(x!.i1-x!.i0)) * y!.i0

closest :: R.Source r Double => Vec r -> Double -> Int
closest x x0 = runIdentity $ do
  -- Differences of xs to x0 with indices
  let dx = R.traverse x id (diffWithIndex (R.extent x) x0)
  fst <$> (!Z) <$> R.foldP minWithIndex (-1,1e100) dx

diffWithIndex :: (Shape sh) => sh -> Double -> (sh -> Double) -> sh -> (Int, Double)
diffWithIndex ext x0 l sh = (R.toIndex ext sh,abs $ l sh - x0)

minWithIndex :: (Int, Double) -> (Int, Double) -> (Int, Double)
minWithIndex (sh1,dx1) (sh2,dx2) = if dx1<dx2 then (sh1,dx1) else (sh2,dx2)

{-
let x = vecFromList [1.0,2.0,3.0,4.0,5.0,7.0,10.0]
let y = vecFromList[5.0,5.0,7.0,5.8,6.0,7.0, 7.5]
let x0=2.8
linearInterpolation x y x0
logLinearInterpolation x y x0
-}

v ^! i = v A.! i

--- |Given coordinates xa[1..n] and ya[1..n], and given a value x, this routine
--- returns a value y, and an error estimate dy. If P(x) is the polynomial of
--- degree n >= 1 such that P(xai) = yai; i=1,...,n, then the returned value
--- is y = P(x) together with an error estimate. The list of points shall have
--- the following form:
---      ps = array (1,3) [(1,(1,2)),(2,(3,4)),(3,(7,1))] :: Array Int (Double,Double)
polynomialInterpolation :: A.Array Int (Double,Double) -> Double -> (Double,Double)
polynomialInterpolation a x = let l=(snd (a^!ic)):step n0 (ic-1)
                              in (sum l,last l)
    where   (n0,n1) = A.bounds a
            cmpByX (_,(x1,_)) (_,(x2,_)) = compare (abs (x-x1)) (abs (x-x2))
            closest a x = L.minimumBy cmpByX (A.assocs a)
            ic = fst $ closest a x  -- index of closest
            c m i = let xi = fst $ a^!i
                        xipm = fst $ a^!(i+m)
                    in if xi==xipm then equalXs
                       else if m==n0 then (xi-x)*((snd $ a^!(i+1)) - (snd $ a^!i))/(xi-xipm)
                       else (xi-x)*((c (m-1) (i+1)) - (d (m-1) i))/(xi-xipm)
            d m i = let xi = fst $ a^!i
                        xipm = fst $ a^!(i+m)
                    in if xi==xipm then equalXs
                       else if m==n0 then (xipm-x)*((snd $ a^!(i+1)) - (snd $ a^!i))/(xi-xipm)
                       else (xipm-x)*((c (m-1) (i+1)) - (d (m-1) i))/(xi-xipm)
            step m ic | m==n1-1 = if ic==0 then [c m 1] else [d m ic]   -- last step: c or d
                      | 2*(ic-n0)<n1-m = (c m (ic+1)):step (m+1) ic     -- c step
                      | otherwise = (d m ic):step (m+1) (ic-1)          -- d step
            equalXs = error "polynomialInterpolation: Points must not have equal x-coordinates."
