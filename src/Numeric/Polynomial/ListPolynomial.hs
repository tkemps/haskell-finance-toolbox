{-# LANGUAGE FlexibleInstances,IncoherentInstances #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  CRP.Polynomial
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

module Numeric.Polynomial.ListPolynomial (
  -- *The Polynomial abstract data type
  Polynomial,
  -- *Constructors
  fromList,zero,one,
  -- *Conversion
  toList,
  -- *Printing
  ppP,ppP',
  -- *Operations
  order,infinite,eval,eval1,takeP,(@!),pscale,psum,pprod,
  -- *Parameters
  nMax
  ) where

import Control.Applicative
import Text.Printf

-- |A @Polynomial a@ represents a polynomial in one variable, e.g. 1+2*x^2+3*x^4,
-- with coefficients of type a and Int exponents >= 0. A polynomial can be constructed
-- with 'fromList' or in special cases with 'one' and 'zero'.
--
-- The implementation of the instances Fractional and Floating create infinite lazy
-- Polynomials that should be handled with appropriate care, e.g. @(order (exp p))@

-- does not return a value. There is no optimization that would recognize that e.g.
-- @(exp one)@ is in fact a finite polynomial.
--
-- Only the following functions from Floating are implemented: '**', 'log', 'exp', 'pi'. The following functions from Floating are not implemented: 'logBase', 'sin', 'cos',
-- 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh'.
data Polynomial a = P [a] Bool deriving Eq

instance (Show a) => Show (Polynomial a) where
  show = shP show

instance Show (Polynomial Double) where
  show = shP shDs

-- |Standard maximum number of items recognized by 'eval1', 'show' and 'ppP'.
nMax :: Int
nMax = 300

more :: Bool -> String
more inf = "..."++(if inf then "(infinite)" else "")

-- |Just like 'ppP' but with an explicit upper limit for the number of terms printed.
ppP' :: Int -> Polynomial Double -> String
ppP' m (P l inf) = (foldl1 (++) (map (uncurry showTerm) v1)) ++ more'
  where v = filter (\(a,_) -> a/=0) (zip l [0..])
        (v1,v2) = splitAt m v
        more' = if null(v2) then "" else (" + "++more inf)
        showTerm a n = " "++(shD a)++" * x^"++(show n)

shP :: ([a] -> String) -> Polynomial a -> String
shP showFun (P v inf) = "Polynomial <<<"++(showFun v1)++more'++">>>"
  where (v1,v2) = splitAt nMax v
        more' = if null(v2) then "" else (", "++more inf)

shDs :: [Double] -> String
shDs xs = foldr1 (\a b -> a++", "++b) (map shD xs)

shD :: Double -> String
shD = printf "%+.4g"

-- |Pretty print a polynomial in the format @a * x^0 + b * x^1 + c * x^2@ etc. up
-- to 'nMax' terms.
ppP :: Polynomial Double -> [Char]
ppP = ppP' nMax

instance Functor Polynomial where
  fmap f (P l inf) = fromList (map f l) inf

-- |Return 'Just' the order of the polynomial, i.e. the highest exponent, or
-- 'Nothing' if it is infinite. The function does not check whether coefficient
-- of the term with the highest exponent is in fact not equal to zero.
order :: Polynomial a -> Maybe Int
order (P l inf) = if inf then Nothing else Just ((length l)-1)

-- |Return whether or not the polynomial is inifinte.
infinite :: Polynomial a -> Bool
infinite (P _ inf) = inf

-- |An infinite lazy list of integer powers of a calculated by successive
-- multiplication with @a@.
powers :: Num a => a -> [a]
powers x = pow 1
  where pow p = p:(pow (p*x))

-- |Evaluate the polynomial @p@ at a certain point @x@.
eval :: Num a =>
        Int            -- ^The maximum order until that evaluation shall continue
        -> Polynomial a -- ^The polynomial p
        -> a            -- ^The point @x@ where to evaluate the polynomial
        -> a
eval nMax p@(P l _) x = sum (map (uncurry (*)) (zip (take nMax l) (take nMax (powers x))))

-- |Evaluation of a polynomial like 'eval' but up to an implicit maximum order 'nMax'.
eval1 :: Num a => Polynomial a -> a -> a
eval1 = eval nMax

-- |Build a polynomial from a list of coefficients. E.g. @fromList [1,0,2]@ represents
-- @1+x^2@.
fromList :: [a]       -- ^List of coefficients (may be a lazy infinite list)
            -> Bool   -- ^Is the list infinite?
            -> Polynomial a
fromList = P

-- |Return the list of coefficients.
toList :: Polynomial t -> [t]
toList (P l _) = l

-- |Equivalent to @fromList [0] False@
zero :: Num a => Polynomial a
zero = pure 0

-- |Equivalent to @fromList [1] False@
one :: Num a => Polynomial a
one = pure 1

-- |Return the first coefficients.
takeP :: Int -> Polynomial a -> Polynomial a
takeP n (P l _) = fromList (take n l) False

build :: Int -> (Int -> a) -> [a]
build n f = map f [0..(n-1)]

(!!!) :: (Num a, Num i,Eq i) => [a] -> i -> a
[] !!! _ = error "Empty list cannot be indexed."
(v0:vs) !!! 0 = v0
(_:[]) !!! _ = 0
(_:vs) !!! n = vs !!! (n-1)

-- |Return a specific coefficient or 0 if the polynomial is shorter than the
-- specified exponent.
(@!) :: (Num a, Num i,Eq i) => Polynomial a -> i -> a
(P l _) @! n = l!!!n

rest :: [a] -> [a] -> [a]
rest v w | lv<lw = drop lv w
         | lv>lw = drop lw v
         | otherwise = []
  where lv = length v
        lw = length w

-- |Lazy sum of polynomials.
psum :: Num a => [Polynomial a] -> Polynomial a
psum ps =  foldl1 (+) ps

-- |Lazy product of polynomials.
pprod :: Num a => [Polynomial a] -> Polynomial a
pprod ps =  foldl1 (*) ps

-- |Multiply each coefficient of the polynomial with a scalar.
pscale :: Num a => a -> Polynomial a -> Polynomial a
pscale x (P l inf) = P (map (*x) l) inf

pneg :: Num a => Polynomial a -> Polynomial a
pneg (P v inf) = P (map negate v) inf

padd :: Num a => Polynomial a -> Polynomial a -> Polynomial a
padd (P v infV) (P w infW) = P ((map (uncurry (+)) (zip v w)) ++ (rest v w))
                               (infV || infW)

pmul :: Num a => Polynomial a -> Polynomial a -> Polynomial a
pmul (P v infV) (P w infW) = P (prod v w) (infV || infW)
  where prod1 v w n = sum [v!!!(n-i)*w!!!i | i<-[0..n]]
        prod v w = build ((length v-1)+(length w-1)) (prod1 v w)

pdiv :: Fractional a => Polynomial a -> Polynomial a -> Polynomial a
pdiv (P v _) (P w _) = fromList div True
  where div1 n = (v!!!n - sum [w!!!i * (div!!(n-i)) | i<-[1..n]]) / (w!!!0)
        div = map div1 [0..]

pdiv' :: Fractional a => Int -> Polynomial a -> Polynomial a -> Polynomial a
pdiv' n p q = takeP n (pdiv p q)

divI :: (Integral a1, Integral a2, Fractional a) => a1 -> a2 -> a
divI i j = (fromIntegral i)/(fromIntegral j)

-- |The function call @pexp p@ creates a polynomial that represents the
-- Tylor expansion of @exp p@. The result is an infinite lazy list.
pexp :: Floating a => Polynomial a -> Polynomial a
pexp p = fromList e True
  where e1 0 = exp (p@!0)
        e1 n = sum [(divI j n)*p@!j*e!!(n-j) | j <- [1..n]]
        e = map e1 [0..]

plog :: Floating a => Polynomial a -> Polynomial a
plog p = fromList l True
  where l1 0 = log (p@!0)
        l1 n = (p@!n - sum [(divI j n) * p@!(n-j) * l!!j | j <- [1..(n-1)]])/p@!0
        l = map l1 [0..]

instance (Num a) => Num (Polynomial a) where
  (+) = padd
  (*) = pmul
  negate = pneg
  abs = fmap abs
  signum = fmap signum
  fromInteger n = fromList [fromInteger n] False

instance (Fractional a) => Fractional (Polynomial a) where
  (/) = pdiv
  fromRational x = fromList [fromRational x] False

instance (Floating a) => Floating (Polynomial a) where
  pi = fromList [pi] False
  exp = pexp
  log = plog
  a ** b = pexp (b * (plog a))
  logBase = undefined
  sin = undefined
  cos = undefined
  tan = undefined
  asin = undefined
  acos = undefined
  atan = undefined
  sinh = undefined
  cosh = undefined
  tanh = undefined
  asinh = undefined
  acosh = undefined
  atanh = undefined

instance Applicative Polynomial where
  pure x = P [x] False
  (P fs inf1) <*> (P as inf2) = P (map (uncurry ($)) (zip fs as)) (inf1 && inf2)
