{-# LANGUAGE ForeignFunctionInterface,NamedFieldPuns,RecordWildCards,ImplicitParams,CPP #-}
module FFI.Tools
       (convD,multiSplitAt) where

import qualified Control.Exception as E
import Data.List
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU

convD :: (Real a,Fractional b) => a -> b
convD = fromRational . toRational

multiSplitAt :: Int -> [a] -> [[a]]
multiSplitAt n l = mat' l []
  where mat' [] acc = acc
        mat' l acc = let (a,b) = splitAt n l
                     in mat' b (acc++[a])

