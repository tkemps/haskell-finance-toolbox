{-# LANGUAGE FlexibleContexts #-}
module Pricing.AmericanOptions (americanOption,americanOption',
                                americanCall,americanPut,
                                callPayoff,putPayoff) where

import           Control.Monad.Identity
import           Data.Array.IArray (listArray)
import           Data.Array.Repa ((!))
import           Data.Array.Repa (DIM1,DIM2)
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import           Data.Array.Repa.Shape
import           Data.List  (sortBy)

type Vec r = R.Array r DIM1 Double
type VecD = Vec R.D
type VecU = Vec R.U

type Mat r = R.Array r DIM2 Double
type MatU = Mat R.U
type MatD = Mat R.D

minVec :: (R.Source r Double) => Vec r -> Vec r -> VecD
minVec v w = R.zipWith min v w

closest :: (Ord a, Num a) => a -> [a] -> [(Int, a)]
closest z l = let l' = zip [0..] (map (\x -> abs (x-z)) l)
                  cmp (_,x) (_,y) = compare x y
              in sortBy cmp l'

americanCall vol d r k t1 n x = americanOption (callPayoff k) vol d r k t1 n x
americanPut vol d r k t1 n x = americanOption (putPayoff k) vol d r k t1 n x

americanOption :: (Double -> Double) -> Double -> Double -> Double -> Double -> Double -> Int -> Double
             -> (Double, Int, Double, Int, Double, VecU)
americanOption payoff vol d r k t1 n x =
  let (n,ds,m,d,prices) = americanOption' payoff vol d r k t1 n
      [(i1,_),(i2,_)] = take 2 $ closest x ([ds*fromIntegral i | i<-[0..n]])
      j1 = min i1 i2
      j2 = max i1 i2
      p1 = prices!(ix1 j1)
      p2 = prices!(ix1 j2)
      s1 = ds * fromIntegral j1
      s2 = ds * fromIntegral j2
      price = (p2-p1)/(s2-s1)*(x-s1) + p1
  in (price,n,ds,m,d,prices)

vectorFromList l = let n = length l in R.fromListUnboxed (ix1 n) l

-- |Price matrix of an American call option with the following parameters:
-- payoff = payoff function with price and strike as parameters
-- vol = Volatility
-- d = dividend yield
-- r = interest rate
-- k = strike
-- t1 = time to expiry
-- n = number of asset value steps
-- example
-- let (n,dS,m,dt,prices) = americanOption' callPayoff 0.2 0.03 0.05 100 3 20
americanOption' :: (Double -> Double) -> Double -> Double -> Double -> Double -> Double -> Int -> (Int,Double,Int,Double,VecU)
americanOption' payoff vol d r k t1 n = runIdentity $ do
  let dS = 2*k/fromIntegral n
      dt' = 0.9/(fromIntegral n)^2/vol^2
      m = (floor (t1/dt')) + 1
      dt = t1/fromIntegral m
      s = R.fromListUnboxed (ix1 (n+1)) [dS * fromIntegral i | i<-[0..n]]
      po = R.map payoff s
  let calc = linkM m $ \pr -> do
        let pr' = timeStep vol d r dS dt n po pr
        pr'' <- R.computeUnboxedP pr'
        return pr''
  po' <- R.computeUnboxedP po
  prices <- calc po'
  return (n,dS,m,dt,prices)

timeStep :: Double -> Double -> Double -> Double -> Double -> Int -> VecD -> VecU -> VecD
timeStep vol d r dS dt n po v =
  let prices = R.traverse v id (previousPrice vol d r dS dt n)
      pricesWithCapAtPayoff = minVec prices po
  in pricesWithCapAtPayoff

previousPrice :: Double -> Double -> Double -> Double -> Double -> Int -> (DIM1 -> Double) -> DIM1 -> Double
previousPrice vol d r dS dt n v i@(Z :. ii) 
  | ii==0 = let x = v i
            in x*(1-r*dt)
  | ii==n = let y1 = v (addDim i (ix1 (-1)))
                y2 = v (addDim i (ix1 (-2)))
            in 2*y1-y2
  | otherwise = let vp = v (addDim i (ix1 1))
                    v0 = v i
                    vm = v (addDim i (ix1 (-1)))
                    s = dS * fromIntegral ii
                    delta = (vp-vm)/2/dS
                    gamma = (vp-2*v0+vm)/(dS^2)
                    theta = -0.5*vol^2 * s^2 * gamma
                            - (r-d) * s * delta + r*v0
                in v0-theta*dt

linkM :: (Monad m) => Int -> (a -> m a) -> (a -> m a)
linkM 0 _ = return
linkM 1 f = f
linkM n f
  | n>1 = \x -> f x >>= linkM (n-1) f
  | otherwise = error $ "Cannot compose "++show n++" times."

callPayoff :: Double -> Double -> Double
callPayoff k x = max (x-k) 0

putPayoff :: Double -> Double -> Double
putPayoff k x = max (k-x) 0

{-
-- example:
-- let (n,dS,m,dt,prices) = americanCall 0.2 0.03 0.05 100 1 20
-- R.toList $ R.slice prices (Any :. All :. (m::Int))
-- R.toList $ R.slice prices (Any :. All :. (0::Int))
americanCall' :: Double -> Double -> Double -> Double -> Double -> Int
                 -> (Int,Double,Int,Double,MatU)
americanCall' vol d r k t1 n = runST $ do
  let dS = 2*k/fromIntegral n
      dt' = 0.9/(fromIntegral n)^(2::Int)/vol^(2::Int)
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
          gamma = (vp-2*v0+vm)/dS^(2::Int)
          theta = -0.5*vol^(2::Int) * s!(Z :. i)^(2::Int) * gamma - (r-d) * s!(Z :. i) * delta + r*v0
      writeArray v (i,j) (v0-theta*dt)
    x <- readArray v (0,j-1)
    writeArray v (0,j) (x*(1-r*dt))
    y1 <- readArray v (n-1,j)
    y2 <- readArray v (n-2,j)
    writeArray v (n,j) (2*y1-y2)
    forM_ [0..n] $ \ i -> do
      vij <- readArray v (i,j)
      writeArray v (i,j) (max vij (payoff!(Z :. i)))
  l <- getElems v
  return (n,dS,m,dt,R.fromListUnboxed (Z :. n+1 :. m+1) l)

americanCall' :: Double -> Double -> Double -> Double -> Double -> Int
                 -> (Int,Double,Int,Double,R.Array R.U ((Z :. Int) :. Int) Double)
americanCall' vol d r k t1 n = runST $ do
  v <- newArray ((0,0),(n,m)) 0 :: ST s (STArray s (Int,Int) Double)
  forM_ [0..n] $ \ i -> writeArray v (i,0) (payoff!(Z :. i))
  forM_ [1..m] (timeStep v)
  l <- getElems v
  return (n,dS,m,dt,R.fromListUnboxed (Z :. n+1 :. m+1) l)

  where dS = 2*k/fromIntegral n
        dt' = 0.9/(fromIntegral n)^(2::Int)/vol^(2::Int)
        m = (floor (t1/dt')) + 1 :: Int
        dt = t1/fromIntegral m
        s = R.fromListUnboxed (Z :. (n+1::Int)) [dS * fromIntegral i | i<-[0..n]]
        payoff = R.map (\x -> max (x-k) 0) s

        timeStep v j = do
          forM_ [1..n-1] (priceStep v j)
          x <- readArray v (0,j-1)
          writeArray v (0,j) (x*(1-r*dt))
          y1 <- readArray v (n-1,j)
          y2 <- readArray v (n-2,j)
          writeArray v (n,j) (2*y1-y2)
          forM_ [0..n] $ \ i -> do
            vij <- readArray v (i,j)
            writeArray v (i,j) (max vij (payoff!(Z :. i)))

        priceStep v j i = do
          vp <- readArray v (i+1,j-1)
          vm <- readArray v (i-1,j-1)
          v0 <- readArray v (i  ,j-1)
          let delta = (vp-vm)/2/dS
              gamma = (vp-2*v0+vm)/dS^(2::Int)
              theta = -0.5*vol^(2::Int) * s!(Z :. i)^(2::Int) * gamma - (r-d) * s!(Z :. i) * delta + r*v0
          writeArray v (i,j) (v0-theta*dt)
-}
