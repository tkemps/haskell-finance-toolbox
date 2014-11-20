module Control.Monad.Util where

linkM :: (Monad m) => Int -> (a -> m a) -> (a -> m a)
linkM 0 _ = return
linkM 1 f = f
linkM n f
  | n>1 = \x -> f x >>= linkM (n-1) f
  | otherwise = error $ "Cannot compose "++show n++" times."
