{-# LANGUAGE FlexibleContexts #-}
module Numeric.Interpolation.Polynomial (
  linearInterpolation, logLinearInterpolation,cubicInterpolation
  ) where

import           Control.Applicative
import           Control.Monad.Identity
import           Control.Monad.ST
import           Data.Array.Repa ((!),DIM1,DIM2)
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import qualified Data.Array.Repa.Operators.Mapping as R
import           Data.Array.Repa.Shape
import           Data.Array.Repa.Slice
import           Data.List  (sortBy)
import qualified Data.Vector.Unboxed as V hiding (length)
import qualified Data.Vector.Unboxed.Mutable as V
import System.IO.Unsafe

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

