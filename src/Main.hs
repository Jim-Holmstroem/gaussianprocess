{-# LANGUAGE FlexibleContexts #-}

module Main where

import Graphics.Plot

import Data.Packed.Vector
import Data.Packed.Matrix
import Numeric.Container
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.HMatrix (tr)

import System.Environment (getArgs)



-- TODO (a->b->c) -> Vector a -> Vector b -> Matrix c
outerOp :: (Element a, Element b) => (Matrix a -> Matrix b -> Matrix c) -> Vector a -> Vector b -> Matrix c
outerOp k x x' = k xp xq
    where xp = repmat (asColumn x) 1 (dim x')
          xq= repmat (asRow x') (dim x) 1



k_SE sigma l x x' = (sigma^2*) $ exp $ (-1/(2*l^2))*dx^2
    where dx = outerOp (-) x x'

k = k_SE 3 0.3



inferenceTrivial k x y ySigma2 domain = (mean, variance)
    where mean = k domain x <> invkxx <> y
          variance = k domain domain - (kdomainx <> invkxx <> k x domain)
          invkxx = inv $ k x x + diag ySigma2
          kdomainx = k domain x

inferenceCholesky k x y ySigma2 domain = (mean, variance)
    where mean = (tr kxdomain) <> alpha
          variance = k domain domain - (tr v <> v)
          alpha = tr l <\> (l <\> y)
          l = chol $ k x x + diag ySigma2
          v = l <\> kxdomain
          kxdomain = k x domain





plotProcessDistribution :: Vector Double -> Vector Double -> Vector Double -> IO ()
plotProcessDistribution domain m v = mplot [domain, m, m + v, m - v, m + 2 * v, m - 2 * v]


main :: IO ()
main = do
    datfile:[] <- getArgs
    dat <- fmap readMatrix $ readFile datfile
    let [x, y] = toColumns dat -- [x,y] or [x,f]
    let domain = linspace 128 (0.5, 10.5)

    let n = 1
    let (x',y') = (vjoin $ replicate n x, vjoin $ take n $ repeat y)

    let (m, v') = inferenceTrivial k x' y' (constant 1 $ dim y') domain
    let v = takeDiag v'

    plotProcessDistribution domain m v
