module Numeric.Root where

import Data.Maybe
import Numeric.Machine

rootRidder :: (Double -> Double) -> Double -> Double -> Maybe Double -> Maybe Double
rootRidder f x1 x2 accuracyGoal = if signum (f x1) == signum (f x2) then Nothing else rootRidder1 (x1,f x1) (x2,f x2) Nothing
    where accGoal = fromMaybe eps accuracyGoal
          rootRidder1 (x1,f1) (x2,f2) acc = if s == 0 then Nothing
                                            else if f4 == 0 then Just x4
                                            else if accuracy<=accGoal then Just x4
                                            else if worseAccuracy accuracy acc then Just x4
                                            else rootRidder1 (x1',f1') (x2',f2') (Just accuracy)
                where x3 = (x1+x2)*0.5
                      f3 = f x3
                      s = sqrt (f3*f3-f2*f1)
                      x4 = x3 +(x3-x1)*(signum (f1-f2))*f3/s
                      f4 = f x4
                      (x1',f1',x2',f2') = if (signum f3)/=(signum f4) then (x3,f3,x4,f4)
                                          else if (signum f1)/=(signum f4) then (x1,f1,x4,f4)
                                          else if (signum f2)/=(signum f4) then (x4,f4,x2,f2)
                                          else error "rootRidder: We should never come here."
                      accuracy = abs(x1'-x2')
                      worseAccuracy thisAcc (Just lastAcc) = thisAcc>=lastAcc
                      worseAccuracy _ Nothing = False
