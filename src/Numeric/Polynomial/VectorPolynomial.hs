{-# LANGUAGE FlexibleInstances,IncoherentInstances,ImplicitParams,RankNTypes #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  CRP.VectorPolynomial
-- Copyright   :  (c) SKS Unternehmensberatung GmbH & Co. KG, 2012
-- License     :  see the file libraries/base/LICENSE
-- 
-- Maintainer  :  christian.haggert@sks-ub.de
-- stability   :  experimental
-- Portability :  portable
--
-- Definition of polynomials and corresponding operations. The polynomials
-- may be infinite since lazy lists are used for the implementation.
--
-----------------------------------------------------------------------------

module Numeric.Polynomial.VectorPolynomial (
  -- *The Polynomial abstract data type
  Polynomial,
  -- *Constructors
  fromList,fromDouble,zero,one,
  -- *Conversion
  toList,toVector,
  -- *Printing 
  ppP,
  -- *Operations
  order,eval,(.!),(.+),(.-),(./),(.*),pexp,plog,pscale,psum,pprod,pmap,pnegate
  ) where

import Prelude hiding (fromInteger)
import qualified Prelude as P (fromInteger)

import Control.Applicative
import Control.Monad.ST
import Data.Maybe
import Data.STRef
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import Data.Vector.Unboxed (Vector,(!),(!?))
import Text.Printf

-- |A @Polynomial a@ represents a polynomial in one variable, e.g. 1+2*x^2+3*x^4, 
-- with coefficients of type a and Int exponents >= 0. A polynomial can be constructed
-- with 'fromList' or in special cases with 'one' and 'zero'.
--
-- The implementation of the instances Fractional and Floating create finite
-- Polynomials according to the parameter maxOrder of the polynomial.
--
-- Only the following functions from Floating are implemented: '**', 'log', 'exp'.
data Polynomial = P {toVector :: Vector Double} deriving (Eq,Show)

-- |Pretty print a polynomial in the format @a * x^0 + b * x^1 + c * x^2@ etc.
ppP :: (?nMax :: Int) => Polynomial -> String
ppP (P v) = foldl1 (++) (map (uncurry showTerm) ts)
  where ts = filter (\(a,_) -> a/=0) (zip (VU.toList v) [0 .. ?nMax])
        showTerm a n = " "++(shD a)++" * x^"++(show n)

shDs :: [Double] -> String
shDs xs = foldr1 (\a b -> a++", "++b) (map shD xs)

shD :: Double -> String
shD = printf "%+.4g"

-- |Apply a function to each coefficient.
pmap :: (Double -> Double) -> Polynomial -> Polynomial
pmap f (P v) = P (VU.map f v)

-- |Return 'Just' the order of the polynomial, i.e. the highest exponent, or 
-- 'Nothing' if it is infinite. The function does not check whether coefficient 
-- of the term with the highest exponent is in fact not equal to zero.
order :: Polynomial -> Int
order (P v) = VU.length v

-- |An infinite lazy list of integer powers of a calculated by successive 
-- multiplication with @a@.
powers :: Double -> [Double]
powers x = pow 1
  where pow p = p:(pow (p*x))

-- |Evaluate the polynomial @p@ at a certain point @x@.
eval :: Polynomial -- ^The polynomial p
        -> Double            -- ^The point @x@ where to evaluate the polynomial
        -> Double
eval p@(P v) x = sum (map (uncurry (*)) (zip (VU.toList v) (powers x)))

-- |Build a polynomial from a list of coefficients. E.g. @fromList [1,0,2]@ represents
-- @1+x^2@.
fromList :: (?nMax :: Int)
            => [Double]       -- ^List of coefficients (may be a lazy infinite list)
            -> Polynomial
fromList l = P (VU.fromList (take ?nMax l))

-- |Return the list of coefficients.
toList :: Polynomial -> [Double]
toList = VU.toList . toVector

fromDouble :: Double -> Polynomial
fromDouble x = P (VU.fromList [x])

-- |Equivalent to @fromList [0]@
zero :: Polynomial
zero = fromDouble 0

-- |Equivalent to @fromList [1]@
one :: Polynomial
one = fromDouble 1

-- |Return a specific coefficient or 0 if the polynomial is shorter than the 
-- specified exponent.
(.!) :: (?nMax::Int) => Polynomial -> Int -> Double
(P v) .! n = if n > ?nMax then 0 else fromMaybe 0 (v!?n)

pnegate :: Polynomial -> Polynomial
pnegate = pmap negate

-- |Multiply each coefficient of the polynomial with a scalar.
pscale :: Double -> Polynomial -> Polynomial
pscale x = pmap (*x)

rest :: [a] -> [a] -> [a]
rest v w | lv<lw = drop lv w
         | lv>lw = drop lw v
         | otherwise = []
  where lv = length v
        lw = length w

(.+) :: Polynomial -> Polynomial -> Polynomial
(P v1) .+ (P v2) = P (VU.fromList $ (map (uncurry (+)) (zip l1 l2)) ++ (rest l1 l2))
  where l1 = VU.toList v1
        l2 = VU.toList v2

p .- q = p .+ (pnegate q)

-- |Lazy sum of polynomials.
psum :: (?nMax::Int) => [Polynomial] -> Polynomial
psum ps =  foldl1 (.+) ps

(!!!) :: Vector Double -> Int -> Double
v !!! n = fromMaybe 0 (v!?n)

build :: Int -> (Int -> Double) -> Vector Double
build n f = VU.fromList $ map f [0..(n-1)]

(.*) :: (?nMax::Int) => Polynomial -> Polynomial -> Polynomial
(P v1) .* (P v2) = P (prod v1 v2)
  where prod1 :: Vector Double -> Vector Double -> Int -> Double
        prod1 v w n = sum [v!!!(n-i)*w!!!i | i<-[0 .. ?nMax]]
        prod :: Vector Double -> Vector Double -> Vector Double
        prod v w = build ((VU.length v-1)+(VU.length w-1)) (prod1 v w)

-- |Lazy product of polynomials.
pprod :: (?nMax::Int) => [Polynomial] -> Polynomial
pprod ps =  foldl1 (.*) ps

-- |Series expansion up to order @?nMax@ of the quotient of the polynomials.
(./) :: (?nMax::Int) => Polynomial -> Polynomial -> Polynomial
(P v) ./ (P w) = fromList div
  where div1 n = (v!!!n - sum [w!!!i * (div!!(n-i)) | i<-[1 .. n]]) / (w!!!0)
        div = map div1 [0 .. ?nMax]

divI :: (Integral a1, Integral a2, Fractional a) => a1 -> a2 -> a
divI i j = (fromIntegral i)/(fromIntegral j)

-- |The function call @pexp p@ creates a polynomial that represents the
-- Tylor expansion of @exp p@. The result is an infinite lazy list.
pexp :: (?nMax::Int) => Polynomial -> Polynomial
pexp p = fromList e
  where e1 0 = exp (p.!0)
        e1 n = sum [(divI j n)*p.!j*e!!(n-j) | j <- [1..n]]
        e = map e1 [0 .. ?nMax]

{-
plog :: (?nMax::Int) => Polynomial -> Polynomial
plog (P v) = P . runST $ plog' v

plog' :: VU.Vector Double -> forall s.ST s (VU.Vector Double)
plog' v = do
  n <- newSTRef 0
  return v
-}
plog :: (?nMax::Int) => Polynomial -> Polynomial
plog p = fromList l
  where l1 0 = log (p.!0)
        l1 n = (p.!n - sum [(divI j n) * p.!(n-j) * l!!j | j <- [1..(n-1)]])/p.!0
        l = map l1 [0 .. ?nMax]

test n = let ?nMax=n 
         in do
            let p1 = fromList [1,2,3]
            let p2 = fromList [1,2,3,4,5]
            putStrLn (show (p1./p2))
            
