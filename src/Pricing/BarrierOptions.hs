module BarrierOptions (
  upAndOutCall,upAndInCall,
  downAndOutCall,downAndInCall,
  downAndOutPut,downAndInPut,
  upAndOutPut,upAndInPut
  ) where

import Numeric.Distribution.Normal

upAndOutCall :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
upAndOutCall s t1 vol r d k b = s*ddt1*(nd1' - nd3' - bb'*(nd6' - nd8'))
                                - k*drt1*(nd2' - nd4' - aa'*(nd5' - nd7'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

upAndInCall :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
upAndInCall s t1 vol r d k b = s*ddt1*(nd3' + bb'*(nd6' - nd8'))
                               - k*drt1*(nd4' + aa'*(nd5' - nd7'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

downAndOutCall :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
downAndOutCall s t1 vol r d k b =
  if k > b
  then s*ddt1*(nd1' - bb'*(1 - nd8')) - k*drt1*(nd2' - aa'*(1 - nd7'))
  else s*ddt1*(nd3' - bb'*(1 - nd6')) - k*drt1*(nd4' - aa'*(1 - nd5'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

downAndInCall :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
downAndInCall s t1 vol r d k b =
  if k > b
  then s*ddt1*bb'*(1 - nd8') - k*drt1*aa'*(1 - nd7')
  else s*ddt1*(nd1' - nd3' + bb'*(1 - nd6')) - k*drt1*(nd2' - nd4' + aa'*(1 - nd5'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

downAndOutPut :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
downAndOutPut s t1 vol r d k b = -s*ddt1*(nd3' - nd1' - bb'*(nd8' - nd6'))
                                 + k*drt1*(nd4' - nd2' - aa'*(nd7' - nd5'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

downAndInPut :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
downAndInPut s t1 vol r d k b = -s*ddt1*(1 - nd3' + bb'*(nd8' - nd6'))
                                + k*drt1*(1 - nd4' + aa'*(nd7' - nd5'))
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

upAndOutPut :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
upAndOutPut s t1 vol r d k b =
  if k > b
  then -s*ddt1*(1 - nd3' - bb'*nd6') + k*drt1*(1 - nd4' - aa'*nd5')
  else -s*ddt1*(1 - nd1' - bb'*nd8') + k*drt1*(1 - nd2' - aa'*nd7')
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

upAndInPut :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
upAndInPut s t1 vol r d k b =
  if k > b
  then -s*ddt1*(nd3' - nd1' + bb'*nd6') + k*drt1*(nd4' - nd2' + aa'*nd5')
  else -s*ddt1*bb'*nd8' + k*drt1*aa'*nd7'
  where ddt1 = exp $ -d*t1
        drt1 = exp $ -r*t1
        nd1' = nd1 s t1 vol r d k b
        nd2' = nd2 s t1 vol r d k b
        nd3' = nd3 s t1 vol r d k b
        nd4' = nd4 s t1 vol r d k b
        nd5' = nd5 s t1 vol r d k b
        nd6' = nd6 s t1 vol r d k b
        nd7' = nd7 s t1 vol r d k b
        nd8' = nd8 s t1 vol r d k b
        aa' = aa s vol r d b
        bb' = bb s vol r d b

aa :: Double -> Double -> Double -> Double -> Double -> Double
aa s vol r d b = (b/s)**(-1 + 2*(r-d) / vol^2)

bb :: Double -> Double -> Double -> Double -> Double -> Double
bb s vol r d b = (b/s)**(1 + 2*(r-d) / vol^2)

d1 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d1 s t1 vol r d k b = ((log $ s/k) + (r-d + 0.5*vol^2)*t1) / (vol*(sqrt t1))

d2 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d2 s t1 vol r d k b = ((log $ s/k) + (r-d - 0.5*vol^2)*t1) / (vol*(sqrt t1))

d3 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d3 s t1 vol r d k b = ((log $ s/b) + (r-d + 0.5*vol^2)*t1) / (vol*(sqrt t1))

d4 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d4 s t1 vol r d k b = ((log $ s/b) + (r-d - 0.5*vol^2)*t1) / (vol*(sqrt t1))

d5 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d5 s t1 vol r d k b = ((log $ s/b) - (r-d - 0.5*vol^2)*t1) / (vol*(sqrt t1))

d6 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d6 s t1 vol r d k b = ((log $ s/b) - (r-d + 0.5*vol^2)*t1) / (vol*(sqrt t1))

d7 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d7 s t1 vol r d k b = ((log $ s*k / b^2) - (r-d - 0.5*vol^2)*t1) / (vol*(sqrt t1))

d8 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
d8 s t1 vol r d k b = ((log $ s*k / b^2) - (r-d + 0.5*vol^2)*t1) / (vol*(sqrt t1))

nd1 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd1 s t1 vol r d k b = cdfNormal $ d1 s t1 vol r d k b

nd2 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd2 s t1 vol r d k b = cdfNormal $ d2 s t1 vol r d k b

nd3 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd3 s t1 vol r d k b = cdfNormal $ d3 s t1 vol r d k b

nd4 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd4 s t1 vol r d k b = cdfNormal $ d4 s t1 vol r d k b

nd5 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd5 s t1 vol r d k b = cdfNormal $ d5 s t1 vol r d k b

nd6 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd6 s t1 vol r d k b = cdfNormal $ d6 s t1 vol r d k b

nd7 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd7 s t1 vol r d k b = cdfNormal $ d7 s t1 vol r d k b

nd8 :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
nd8 s t1 vol r d k b = cdfNormal $ d8 s t1 vol r d k b

