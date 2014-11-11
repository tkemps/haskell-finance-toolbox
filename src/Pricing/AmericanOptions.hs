{-# LANGUAGE FlexibleContexts,TypeOperators #-}
module Pricing.AmericanOptions where

import           Control.Exception (throw,ArithException(LossOfPrecision))
import           Control.Monad (forM_,when)
import           Control.Monad.Cont (ContT,runContT,callCC,lift)
import           Control.Monad.ST (ST,runST)
import           Data.Array (Array)
import           Data.Array.IArray (IArray)
import           Data.Array.ST (STArray,MArray,newArray,newListArray,writeArray,
                                readArray,getElems,getBounds,getElems,freeze)
import           Data.Array.Repa ((!))
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import           Data.Array.Repa.Shape

-- |Price matrix of an American call option with the following parameters:
-- vol = Volatility
-- d = dividend yield
-- r = interest rate
-- k = strike
-- t1 = time to expiry
-- n = number of asset value steps
-- example let v = americanCall 0.2 0.03 0.05 100 1 100 :: UArray (Int,Int) Double
americanCall :: Double -> Double -> Double -> Double -> Double -> Int
                -> (Int,Int,R.Array R.U ((Z :. Int) :. Int) Double)
americanCall vol d r k t1 n = runST $ do
  let dS = 2*k/fromIntegral n
      dt' = 0.9/(fromIntegral n)^2/vol^2
      m = (floor (t1/dt')) + 1 :: Int
      dt = t1/fromIntegral m
      s = R.fromListUnboxed (Z :. (n+1::Int)) [dS * fromIntegral i | i<-[0..n]]
      payoff = R.map (\x -> max (x-k) 0) s
  v <- newArray ((0,0),(n,m)) 0 :: ST s (STArray s (Int,Int) Double)
  forM_ [0..n] $ \ i -> writeArray v (i,0) (payoff!(Z :. i))
  forM_ [1..m] $ \ j -> do
    forM_ [1..n-1] $ \ i -> do
      vp <- readArray v (i+1,j-1)
      vm <- readArray v (i-1,j-1)
      v0 <- readArray v (i  ,j-1)
      let delta = (vp-vm)/2/dS
          gamma = (vp-2*v0+vm)/dS^2
          theta = -0.5*vol^2 * s!(Z :. i)^2 * gamma - (r-d) * s!(Z :. i) * delta + r*v0
      writeArray v (i,j-1) (v0-theta*dt)
    x <- readArray v (0,j-1)
    writeArray v (0,j) (x*(1-r*dt))
    y1 <- readArray v (n-1,j)
    y2 <- readArray v (n-2,j)
    writeArray v (n,j) (2*y1-y2)
    forM_ [0..n] $ \ i -> do
      vij <- readArray v (i,j)
      writeArray v (i,j) (max vij (payoff!(Z :. i)))
  l <- getElems v
  return (n,m,R.fromListUnboxed (Z :. n+1 :. m+1) l)
