{-# LANGUAGE FlexibleContexts #-}
module Numeric.Interpolation.Spline (
  cubicSpline
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

cubicSpline :: (R.Source r Double) => Vec r -> Vec r -> Double -> Double
cubicSpline x y = let dx = vecDiff x   -- length n-1
                      dy = vecDiff y   -- length n-1
                      q = R.traverse2 dx dy   -- length n-2
                          (\shx shy -> if shx/=shy
                                       then error "Shapes of x and y do not match"
                                       else addDim shx (ix1 (-1)))
                          (\dx dy i -> let ip1 = addDim i (ix1 1)
                                       in 3 * (dy ip1/dx ip1 - dy i/dx i))
                      (a,b,c) = runST $ do
                        beta <- V.new n :: ST s (V.MVector s Double)
                        delta <- V.new n :: ST s (V.MVector s Double)
                        V.write beta 0 0.0
                        V.write delta 0 0.0
                        forM_ [1..n-2] $ \i -> do
                          beta'im1 <- V.read beta (i-1)
                          delta'im1 <- V.read delta (i-1)
                          let alpha'i = 2 * (dx!.(i-1) + dx!.i) - dx!.(i-1) * beta'im1
                          V.write beta i $ dx!.i / alpha'i
                          V.write delta i $ (q!.(i-1) - dx!.(i-1)*delta'im1) / alpha'i
                        a <- V.new (n-1) :: ST s (V.MVector s Double)
                        b <- V.new n :: ST s (V.MVector s Double)
                        c <- V.new (n-1) :: ST s (V.MVector s Double)
                        V.write b (n-1) 0
                        forM_ (reverse [0..n-2]) $ \i -> do
                          beta'i <- V.read beta i
                          delta'i <- V.read delta i
                          b'ip1 <- V.read b (i+1)
                          let b'i = delta'i - beta'i*b'ip1
                          V.write b i b'i
                          V.write a i $ dy!.i/dx!.i - dx!.i/3*(b'ip1 + 2*b'i)
                          V.write c i $ (b'ip1 - b'i) / (3 * dx!.i)
                        a' <- V.freeze a
                        b' <- V.freeze b
                        c' <- V.freeze c
                        return (R.fromUnboxed (ix1 (n-1)) a',
                                R.fromUnboxed (ix1 n) b',
                                R.fromUnboxed (ix1 (n-1)) c')
                      in evalSpline x y a b c
                    
  where n = R.size (R.extent x)
        vecDiff x = R.traverse x
                       (\sh -> addDim sh (ix1 (-1)))
                       (\x i -> x (addDim i (ix1 1)) - x i)

evalSpline :: R.Source r Double => Vec r -> Vec r -> VecU -> VecU -> VecU -> Double -> Double
evalSpline x y a b c x0 =
  let n = R.size (R.extent x)
      xs = filter (\(_,x)->x<x0) (zip [0..] (R.toList x))  -- to be optimized!
      i = fst $ last $ if null xs then error $ "evalSpline: x="++show x0++" out of range." else xs
      i' = ix1 (if i==n-1 then error  $ "evalSpline: x="++show x0++" out of range." else i)
      dx = x0-x!i'
      dx2 = dx*dx
      dx3 = dx2*dx
  in y!i' + a!i'*dx + b!i'*dx2 + c!i'*dx3

indexed v = R.traverse v id (\x idx@(Z :. i) -> (x idx, i))

{-
let x = vecFromList [1.0,2.0,3.0,4.0,5.0,7.0,10.0]
let y = vecFromList[5.0,5.0,7.0,5.8,6.0,7.0, 7.5]
let f=cubicSpline x y
-}
