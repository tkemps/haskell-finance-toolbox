{-# LANGUAGE FlexibleContexts,FlexibleInstances,UndecidableInstances #-}
module Data.Vector.AddOns where

import Control.Monad
import Control.Monad.Primitive
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import qualified Data.Vector.Generic as VG
import Data.Vector.Generic ((!))

modifyMVector
  :: (PrimMonad m, VUM.Unbox a) =>
     VUM.MVector (PrimState m) a -> Int -> (a -> a) -> m ()
modifyMVector v i f = do
  x <- VUM.read v i
  let y = f x
  VUM.write v i y

binOpMVector
  :: (PrimMonad m, VUM.Unbox a, VUM.Unbox t) =>
     (a -> t -> a)
     -> VUM.MVector (PrimState m) a
     -> VUM.MVector (PrimState m) t
     -> m ()
binOpMVector op v w = do
  let n = VUM.length v
  when (n/=VUM.length w) $ do
    error "binOpMVector: Vectors have different lengths"
  forM_ [0..(n-1)] (\i -> do
    x <- VUM.read v i
    y <- VUM.read w i
    VUM.write v i (x `op` y))

binOpMVector'
  :: (PrimMonad m, VUM.Unbox a, VUM.Unbox a1) =>
     (a -> a1 -> a)
     -> VUM.MVector (PrimState m) a -> VU.Vector a1 -> m ()
binOpMVector' op v w = do
  let n = VUM.length v
  when (n/=VU.length w) $ do
    error "binOpMVector': Vectors have different lengths"
  forM_ [0..(n-1)] $ (\i -> do
    x <- VUM.read v i
    VUM.write v i (x `op` (w!i)))

(+=) :: (PrimMonad m, VUM.Unbox Double) =>
        VUM.MVector (PrimState m) Double -> VUM.MVector (PrimState m) Double
     -> m ()
(+=) = binOpMVector (+)

(-=) :: (PrimMonad m, VUM.Unbox Double) =>
        VUM.MVector (PrimState m) Double -> VUM.MVector (PrimState m) Double
     -> m ()
(-=) = binOpMVector (-)

(*=) :: (PrimMonad m, VUM.Unbox Double) =>
        VUM.MVector (PrimState m) Double -> VUM.MVector (PrimState m) Double
     -> m ()
(*=) = binOpMVector (*)

(\=) :: (PrimMonad m, VUM.Unbox Double) =>
        VUM.MVector (PrimState m) Double -> VUM.MVector (PrimState m) Double
     -> m ()
(\=) = binOpMVector (/)

(+=.) :: (PrimMonad m, VUM.Unbox Double) =>
         VUM.MVector (PrimState m) Double
         -> VU.Vector Double
         -> m ()
(+=.) = binOpMVector' (-)

(-=.) :: (PrimMonad m, VUM.Unbox Double) =>
         VUM.MVector (PrimState m) Double
         -> VU.Vector Double
         -> m ()
(-=.) = binOpMVector' (-)

(*=.) :: (PrimMonad m, VUM.Unbox Double) =>
         VUM.MVector (PrimState m) Double
         -> VU.Vector Double
         -> m ()
(*=.) = binOpMVector' (*)

(\=.) :: (PrimMonad m, VUM.Unbox Double) =>
         VUM.MVector (PrimState m) Double
         -> VU.Vector Double
         -> m ()
(\=.) = binOpMVector' (/)

mapMVector :: (PrimMonad m, VUM.Unbox a) =>
              VUM.MVector (PrimState m) a -> (a -> a) -> m ()
mapMVector v f = do
  let n = VUM.length v
  forM_ [0..(n-1)] (\i -> do
    x <- VUM.read v i
    let y = f x
    VUM.write v i y)

unBox :: VUM.Unbox a => V.Vector a -> VU.Vector a
unBox = VU.fromList . V.toList

vdot :: (VG.Vector v Double) => v Double -> v Double -> Double
vdot v w = VG.sum (VG.zipWith (*) v w)

instance (VG.Vector v Double) => Num (v Double) where
  (+) = VG.zipWith (+)
  (-) = VG.zipWith (-)
  (*) = VG.zipWith (*)
  abs = VG.map abs
  negate = VG.map negate
  signum = VG.map signum
  fromInteger = undefined

