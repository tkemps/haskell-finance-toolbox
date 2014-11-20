{-# LANGUAGE FlexibleContexts #-}
module Pricing.Util (
  Pricing(..),OptionType(..),OptionPricingFunction,genericSensi,
  callPayoff,putPayoff,
  print'
  ) where

import           Control.Monad.Identity
import           Data.Array.IArray (listArray)
import           Data.Array.Repa ((!))
import           Data.Array.Repa (DIM1,DIM2)
import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index
import           Data.Array.Repa.Shape
import           Data.List  (sortBy)
import           Numeric.Distribution.Normal
import           Pricing.BlackScholes
import           System.IO.Unsafe

print' :: Show a => String -> a -> b -> b
print' lbl x a = unsafePerformIO $ do
  putStr (show x)
  putStr "="
  putStrLn lbl
  return a

data OptionType = Call | Put deriving (Eq,Show)

data Pricing = Pricing {
  price :: Double,
  delta :: Double,
  elasticity :: Double,
  theta :: Double,
  rho :: Double,
  gamma :: Double,
  dGammaDVol :: Double,
  gammaP :: Double,
  timeGamma :: Double
  } deriving (Eq,Show)

type OptionPricingFunction = Double        -- ^Price of underlying
                             -> Double        -- ^Strike
                             -> Double        -- ^Remaining maturity in years
                             -> Double        -- ^interest rate
                             -> Double        -- ^cost of carry
                             -> Double        -- ^Volatility
                             -> Double

genericSensi
  :: Double        -- ^Shift of underlying price s
     -> Double        -- ^Shift of vola v
     -> Double        -- ^Shift of time t1
     -> Double        -- ^Shift of interest rate r and cost of carry b
     -> OptionPricingFunction
     -> Double        -- ^Price of underlying
     -> Double        -- ^Strike
     -> Double        -- ^Remaining maturity in years
     -> Double        -- ^interest rate
     -> Double        -- ^cost of carry
     -> Double        -- ^Volatility
     -> Pricing
genericSensi ds dv dt dr f s x t1 r b v = Pricing {
    price = p0,
    delta = (pp - pm)/(2*ds),
    elasticity = (pp - pm)/(2*ds)*s/p0,
    theta = if t1<1/365
            then (priceT 0.00001) - p0
            else (priceT (t1-dt)) - p0,
    rho = (priceRB (r+dr) (b+dr) - priceRB (r-dr) (b-dr))/2,
    gamma = (pp - 2*p0 + pm)/ds^2,
    dGammaDVol = ((priceSV (s+ds) (v+dv)) - 2*(priceSV s (v+dv)) + (priceSV (s-ds) (v+dv))
                   - (priceSV (s+ds) (v-dv)) + 2*(priceSV s (v-dv)) - (priceSV (s-ds) (v-dv)))
                  / (2*dv*ds^2) / 100,
    gammaP = s/100*(pp - 2*p0 + pm)/ds^2,
    timeGamma = (priceT (t1+dt) - 2*p0 + priceT (t1-dt))/dt^2
    }
  where priceS s = f s x t1 r b v
        pp = priceS (s+ds)
        p0 = priceS s
        pm = priceS (s-ds)
        priceSV s v = f s x t1 r b v
        priceT t1 = f s x t1 r b v
        priceRB r b = f s x t1 r b v

callPayoff :: Double -> Double -> Double
callPayoff k x = max (x-k) 0

putPayoff :: Double -> Double -> Double
putPayoff k x = max (k-x) 0
