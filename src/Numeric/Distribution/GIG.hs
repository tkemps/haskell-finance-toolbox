{-# LANGUAGE ScopedTypeVariables,FlexibleContexts,StandaloneDeriving,UndecidableInstances #-}
module Numeric.Distribution.GIG where

import Data.Bits
import Control.Monad
import Data.Maybe
import Numeric.GSL.Special.Bessel
import Numeric.GSL.Integration

-- |g > 0, d > 0, l may in general be real but only l>0 supported by GSL, 
-- support of pdf is x>0
pdfGIG l d g x = if x<=0 then 0 
                 else n * x**(l-1) * (exp (-0.5*(d*d/x+g*g*x)))
  where n = 0.5*(g/d)**l / (bessel_Knu (abs l) (d*g))

meanGIG l d g = d * (bessel_Knu (1+abs l) (d*g)) / (g * (bessel_Knu (abs l) (d*g)))

varGIG l d g = d*d/(g*g) * (b2/b0 - (b1/b0)**2)
               where b k = bessel_Knu (k+abs l) (d*g)
                     b0 = b 0
                     b1 = b 1
                     b2 = b 2

modeGIG l d g = (l-1 + (sqrt((l-1)**2+sqrt(d*g))))/(sqrt g)

cdfGIG l d g x = fst (cdfGIG_e x l d g 1e-12)

cdfGIG_e l d g eps x = if x<0
                       then (0,0)
                       else if x<0.01
                            then (0,3e-26)
                            else integrateQNG eps (pdfGIG l d g) 0 x
