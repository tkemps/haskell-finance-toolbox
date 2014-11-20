module Pricing.BlackScholes (
  callOption,callOptionGBS,callDelta,callGamma,callTheta,callSpeed,callVega,callRho,callRhoD,
  putOption,putOptionGBS,putDelta,putGamma,putTheta,putSpeed,putVega,putRho,putRhoD,
  d1,d2) where

import Numeric.Distribution.Normal

callOptionGBS x vol b r k t1 = callOption x vol (b-r) r k t1
putOptionGBS x vol b r k t1 = putOption x vol (b-r) r k t1

-- |Black-Scholes price of an equity call option.
-- x = asset price
-- vol = volatility of asset price
-- dividend yield
-- r = interest rate
-- k = strike price
-- t1 = expiry
callOption x vol d r k t1 = x*(exp (-d*t1))*(n1p x vol d r k t1)
                            - k*( exp(-r*t1))*(n2p x vol d r k t1)

putOption x vol d r k t1 = -x*(exp(-d*t1))*(n1m  x vol d r k t1)
                           +k*(exp (-r*t1))*(n2m x vol d r k t1)

callDelta x vol d r k t1 = (exp (-d*t1)) * (n1p x vol d r k t1)
putDelta x vol d r k t1 = (exp (-d*t1)) * (-1 + n1p x vol d r k t1)

callGamma x vol d r k t1 = (exp (-d*t1)) * (nxd1 x vol d r k t1) / vol / x / (sqrt t1)
putGamma = callGamma

callTheta x vol d r k t1 = a+b
  where a = -vol*x*(exp (-d*t1)) * (nxd1 x vol d r k t1) / 2 / (sqrt t1)
        b = d*x*(n1p x vol d r k t1)*exp(-d*t1)-r*k*(exp (-r*t1))*(n2p x vol d r k t1)
putTheta x vol d r k t1 = a-b
  where a = -vol*x*(exp (-d*t1)) * (nxd1 x vol d r k t1) / 2 / (sqrt t1)
        b = d*x*(n1m x vol d r k t1)*exp(-d*t1)-r*k*(exp (-r*t1))*(n2m x vol d r k t1)

callSpeed x vol d r k t1 =
  -(exp (-d*t1))*(nxd1 x vol d r k t1)/(vol*vol)/(x*x)/t1*(d1 x vol d r k t1)+vol*(sqrt t1)

putSpeed = callSpeed

callVega x vol d r k t1 = x*(sqrt t1) * (exp (-d*t1)) * (nxd1 x vol d r k t1)
putVega = callVega

callRho x vol d r k t1 = k*t1*(exp (-r*t1))*(n2p x vol d r k t1)

putRho x vol d r k t1 = -k*t1*(exp (-r*t1))*(n2m x vol d r k t1)

callRhoD x vol d r k t1 = -t1*x*(exp (-d*t1))*(n1p x vol d r k t1)

putRhoD x vol d r k t1 = t1*x*(exp (-d*t1))*(n1m x vol d r k t1)

n1p = n 1 d1
n2p = n 1 d2
n1m = n (-1) d1
n2m = n (-1) d2

n s dx x vol d r k t1 = cdfNormal (s*(dx x vol d r k t1))

d1 x vol d r k t1 = ((log (x/k)) + (r-d+0.5*vol*vol)*t1)/vol/(sqrt t1)

d2 x vol d r k t1 = (d1 x vol d r k t1)-vol*(sqrt t1)

nxd1 = nx d1
nxd2 = nx d2
nx dx x vol d r k t1 = 1/(sqrt (2*pi))*(exp (-0.5*(dx x vol d r k t1)*(dx x vol d r k t1)))
