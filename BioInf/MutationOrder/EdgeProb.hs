
module BioInf.MutationOrder.EdgeProb where

import           Control.Arrow (second)
import           Control.Monad (forM_)
import           Data.List (nub,sort)
import           Data.Text (Text,unpack)
import           Data.Vector.Unboxed (Unbox)
import           Numeric.Log
import qualified Data.Text as T
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import           Text.Printf

import           ADP.Fusion.Core
import           ADP.Fusion.EdgeBoundary
import           ADP.Fusion.Set1
import           Data.PrimitiveArray hiding (toList)
import           Diagrams.TwoD.ProbabilityGrid
import           FormalLanguage
import           ShortestPath.SHP.EdgeProb

import           BioInf.MutationOrder.RNA



-- | Before using @aInside@ the @ScoreMat@ needs to be scaled
--
-- TODO the @Edge@ needs to be an @EdgeWithActive@ to get the active bits
-- on the left in the set.

aInside :: Monad m => Landscape -> Temperatur -> SigEdgeProb m (Log Double) (Log Double) (BitSet:.From:.To) Int
aInside Landscape{..} temperature = SigEdgeProb
  { edge = \x e -> x * error "landscape @ e scaled by temperature"
  , mpty = \() -> 1
  , node = \n -> 1
  , fini = \l e f -> l * error "landscape at e scaled by temperature" * f
  , h    = SM.foldl' (+) 0
  }
{-# Inline aInside #-}



type TF1 x = TwITbl Id Unboxed EmptyOk (BS1 First I)      x
type TL1 x = TwITbl Id Unboxed EmptyOk (BS1 Last  I)      x
type EB  x = TwITbl Id Unboxed EmptyOk (EdgeBoundary I)   x

type BF1 x b = TwITblBt Unboxed EmptyOk (BS1 First I)    x Id Id b
type BL1 x b = TwITblBt Unboxed EmptyOk (BS1 Last  I)    x Id Id b
type BEB x b = TwITblBt Unboxed EmptyOk (EdgeBoundary I) x Id Id b



-- | Extract the individual partition scores.

edgeProbPartFun :: Landscape -> Double -> [(EdgeBoundary I, Log Double)]
edgeProbPartFun landscape temperature =
  let n       = mutationCount landscape
      (Z:.sF:.sL:.sZ) = mutateTablesST $ gEdgeProb (aInside landscape temperature)
                          (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 1 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 2 0 EmptyOk (fromAssocs (0 :-> 0)    (0 :-> (n-1))                  0 []))
                          EdgeWithBitSet
                          Singleton
                        :: Z:.TF1 (Log Double):.TL1 (Log Double):.EB (Log Double)
      TW (ITbl _ _ _ pf) _ = sZ
      bs' = assocs pf
      pssum = (Numeric.Log.sum $ Prelude.map snd bs') / (fromIntegral n - 1)
      bs = Prelude.map (second (/pssum)) bs'
  in bs
{-# NoInline edgeProbPartFun #-}

-- | Turn the edge probabilities into a score matrix.

edgeProbScoreMat :: (Unbox t) => ScoreMat t -> [(EdgeBoundary I, Log Double)] -> ScoreMat (Log Double)
edgeProbScoreMat (ScoreMat mat ns) xs' = ScoreMat m ns
  where m = fromAssocs l h 0 xs
        (l,h) = bounds mat
        xs = [ ((Z:.f:.t),p) | (f :-> t, p) <- xs' ]

