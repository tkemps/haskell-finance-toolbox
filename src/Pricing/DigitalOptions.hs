module Pricing.DigitalOptions (
  cashOrNothingCall,cashOrNothingPut
                              ) where

import Numeric.Distribution.Normal

-- |This implementation covers various products as follows: b=r for options on non-dividend paying stock, b=r-q for options on stock or index paying a dividend yield of q, b=0 for options on futures, b=r-rf for currency options (where rf is the rate in the second currency).
cashOrNothingCall ::
  Double     -- ^Asset price (s)
  -> Double   -- ^Strike price (x)
  -> Double   -- ^Cash (k)
  -> Double   -- ^Time to maturity (t1)
  -> Double   -- ^Risk-free rate (r)
  -> Double   -- ^Cost of carry (b)
  -> Double   -- ^Volatility (v)
  -> Double
cashOrNothingCall s x k t1 r b v =
  let d = ((log $ s/x) + (b - 0.5*v^2) * t1) / (v * sqrt t1)
  in k * (exp $ -r*t1) * (cdfNormal d)

-- |This implementation covers various products as follows: b=r for options on non-dividend paying stock, b=r-q for options on stock or index paying a dividend yield of q, b=0 for options on futures, b=r-rf for currency options (where rf is the rate in the second currency).
cashOrNothingPut ::
  Double     -- ^Asset price (s)
  -> Double   -- ^Strike price (x)
  -> Double   -- ^Cash (k)
  -> Double   -- ^Time to maturity (t1)
  -> Double   -- ^Risk-free rate (r)
  -> Double   -- ^Cost of carry (b)
  -> Double   -- ^Volatility (v)
  -> Double
cashOrNothingPut s x k t1 r b v =
  let d = ((log $ s/x) + (b - 0.5*v^2) * t1) / (v * sqrt t1)
  in k * (exp $ -r*t1) * (cdfNormal (-d))
