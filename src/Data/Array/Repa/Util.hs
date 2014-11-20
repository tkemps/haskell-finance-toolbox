module Data.Array.Repa.Util where

import qualified Data.Array.Repa as R
import           Data.Array.Repa.Index

vectorFromList l = let n = length l in R.fromListUnboxed (ix1 n) l

