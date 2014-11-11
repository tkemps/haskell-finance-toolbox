{-# LANGUAGE ForeignFunctionInterface #-}
module Numeric.Knapsack.SubSum (decomp) where

import System.IO.Unsafe (unsafePerformIO)
import Foreign.Marshal.Array (withArray,mallocArray,peekArray)
import Foreign.Ptr (Ptr)

-- int decomp(int n, int *w, int *x, int c)
foreign import ccall "decomp" decomp' :: Int -> Ptr Int -> Ptr Int -> Int -> IO Int

-- |This function takes two arguments w and c and solves a subset-sum problem:
-- @
-- @           maximize   \sum_{j=1}^{n} w_{j} x_{j}
-- @           subject to \sum_{j=1}^{n} w_{j} x_{j} \leq c
-- @                      x_{j} \in \{0,1\}, j = 1,\ldots,n
-- @
-- It returns the solution vector and the optimal objective value
decomp :: [Int] -> Int -> ([Bool],Int)
decomp w c =  unsafePerformIO $ do
                withArray w (\w' -> do
                               x' <- mallocArray n
                               z <- decomp' n w' x' c
                               x <- peekArray n x'
                               return (map (==1) x,z))
    where n = length w
