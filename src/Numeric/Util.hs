{-# LANGUAGE RankNTypes,FlexibleContexts #-}
module Numeric.Util where

import Control.Monad
import Control.Monad.ST
import Data.Array.IArray
import Data.Array.MArray
import Data.Array.ST
import Data.STRef
import qualified Data.Vector.Generic as VG
import Text.Printf
import Statistics.Sample
import Statistics.Quantile

putVec :: forall a e . (IArray a e,Show e) => a Int e -> IO ()
putVec v = forM_ [l..h] (\i->do
                              putStr ((show i)++" ")
                              putStrLn (show (v!i)))
  where (l,h)=bounds v

showMArray :: (Show e,MArray a e m,Ix i) => a i e -> m [String]
showMArray v = do
  xs <- getElems v
  mapM (\x -> return (show x)) xs

binarySearchIArray :: forall a i e.(Ix i,Integral i,Ord e,IArray a e)
                     => a i e -> e -> Maybe i
binarySearchIArray v x = bs i0 i1
  where (i0,i1) = bounds v
        bs l h = if h<l then Nothing
                 else let m = l+(h-l) `div` 2
                      in if v!m>x then bs l (m-1)
                         else if v!m<x then bs (m+1) h
                              else Just m

-- |The array v with n elements constitutes @n-1@ disjoint, consecutives intervals:
-- [v!0;v!1], (v!1;v!2],..., [v!n-2;v!(n-1)] The first interval is closed the
-- others are open on the left and closed on the right side.
-- The function @binaryIntervalSearchIArray@ returns @Just@ the index i such that x 
-- lies within the interval (v!i;v!i+1] for 1<i<n resp. [v!1;v!2]. If
-- no such interval exists since x < v!1 or x > v!n-1 the function returns
-- @Nothing@.
binaryIntervalSearchIArray :: (Ix i,Integral i,Ord e,IArray a e) 
                             => a i e -> e -> Maybe i
binaryIntervalSearchIArray v x = bs i0 i1
  where (i0,i1) = bounds v
        bs l h = if h<l then Nothing
                 else let m = l+(h-l) `div` 2
                      in if v!m>x then bs l (m-1)
                         else if m<i1 && v!(m+1)<x then bs (m+1) h
                              else if m>=i1 then Nothing else Just m

solveTriDiag :: (IArray a Double) =>
     a Int Double -> a Int Double -> a Int Double -> a Int Double -> a Int Double
solveTriDiag a b c r =
  if b!0==0 then error "solveTriDiag: b!0==0, rewrite equations with first equation eliminated."
  else runST $ (do
    u <- newArray (0,n-1) 0 :: ST s (STArray s Int Double)
    writeArray u 0 (r!0/b!0)
    bet <- newSTRef (b!0) :: ST s (STRef s Double)
    gam <- newArray (0,n-1) 0 :: ST s (STArray s Int Double)
    forM_ [1..n-1] 
      (\j->do
          bet' <- readSTRef bet
          writeArray gam j (c!(j-1)/bet')
          gamj <- readArray gam j
          writeSTRef bet (b!j-a!j*gamj)
          bet'' <- readSTRef bet
          when (bet''==0) (error ("solveTriDiag: Zero pivot at column "++(show j)))
          ujm1 <- readArray u (j-1)
          writeArray u j ((r!j-a!j*ujm1)/bet''))
    forM_ [n-2..0] 
      (\j->do
          uj <- readArray u j
          ujp1 <- readArray u (j+1)
          gamjp1 <- readArray gam (j+1)
          writeArray u j (uj-gamjp1*ujp1))
    freeze u)
  where n = snd (bounds a)

analyseSample :: VG.Vector v Double => v Double -> IO ()
analyseSample xs = do
  putStrLn $ printf "Length:          %15d" (VG.length xs)
  putStrLn $ printf "Mean:            %15.2f" (mean xs)
  putStrLn $ printf "StdDev:          %15.2f" (stdDev xs)
  putStrLn $ printf "Skewness:        %15.2f" (skewness xs)
  putStrLn $ printf "Kurtosis:        %15.2f" (kurtosis xs)
  forM_ [0.95::Double,0.98,0.99,0.999] 
    (\q -> do
        putStrLn $ printf "%2.2f%% Qauntile: %15.2f" (100*q) (quantile q))
  where quantile q = continuousBy spss (round (10000*q)) 10000 xs
