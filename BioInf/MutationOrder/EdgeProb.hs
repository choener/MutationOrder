
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
import qualified Data.Vector as V
import qualified Data.HashMap.Strict as HM
import           Data.Bits
import           Debug.Trace

import           ADP.Fusion.Core
import           ADP.Fusion.EdgeBoundary
import           ADP.Fusion.Set1
import           Data.PrimitiveArray hiding (toList)
import           Data.PrimitiveArray.ScoreMatrix
import           Diagrams.TwoD.ProbabilityGrid
import           FormalLanguage
import           ShortestPath.SHP.EdgeProb
import           Data.Vector.Generic.Unstream

import           BioInf.MutationOrder.RNA
import           BioInf.MutationOrder.MinDist (ScaleFunction(..),scaleFunction)



-- | Before using @aInside@ the @ScoreMat@ needs to be scaled
--
-- TODO the @Edge@ needs to be an @EdgeWithActive@ to get the active bits
-- on the left in the set.

aInside :: Monad m => ScaleFunction -> Landscape -> Double -> SigEdgeProb m (Log Double) (Log Double) (Int:.From:.To) Int
aInside scaled Landscape{..} temperature = SigEdgeProb
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                          fene = centroidEnergy frna
                                          tene = centroidEnergy trna
                                      in  traceShow ('E',BitSet fset,f,t,x,fene,tene) $ x * (Exp . negate $ scaleFunction scaled (tene - fene) / s)
  , mpty = \() -> 1
  , node = \n -> let frna = rnas HM.! (BitSet 0)
                     trna = rnas HM.! (BitSet 0 `setBit` n)
                     fene = centroidEnergy frna
                     tene = centroidEnergy trna
                 in  traceShow ('N',n,fene,tene) $ Exp . negate $ scaleFunction scaled (tene - fene) / s
  , fini = \l (fset:.From f:.To t) r -> let frna = rnas HM.! (BitSet fset)
                                            trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                            fene = centroidEnergy frna
                                            tene = centroidEnergy trna
                                        in  traceShow ('F',BitSet fset,f,t,l,r,fene,tene) $ l * r * (Exp . negate $ scaleFunction scaled (tene - fene) / s)
  , h    = SM.foldl' (+) 0
--  , h    = \s -> do v :: V.Vector (Log Double) <- streamToVectorM s
--                    return $ Numeric.Log.sum v
  } where !s = temperature * n * (n-1)
          !n = fromIntegral mutationCount
{-# Inline aInside #-}



type TF1 x = TwITbl Id Unboxed EmptyOk (BS1 First I)      x
type TL1 x = TwITbl Id Unboxed EmptyOk (BS1 Last  I)      x
type EB  x = TwITbl Id Unboxed EmptyOk (EdgeBoundary I)   x

type BF1 x b = TwITblBt Unboxed EmptyOk (BS1 First I)    x Id Id b
type BL1 x b = TwITblBt Unboxed EmptyOk (BS1 Last  I)    x Id Id b
type BEB x b = TwITblBt Unboxed EmptyOk (EdgeBoundary I) x Id Id b



-- | Extract the individual partition scores.

edgeProbPartFun :: ScaleFunction -> Double -> Landscape -> [(EdgeBoundary I, Log Double)]
edgeProbPartFun scaled temperature landscape =
  let n       = mutationCount landscape
      (Z:.sF:.sL:.sZ) = mutateTablesST $ gEdgeProb (aInside scaled landscape temperature)
                          (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 1 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 2 0 EmptyOk (fromAssocs (0 :-> 0)    (0 :-> (n-1))                  0 []))
                          EdgeWithSet
                          Singleton
                        :: Z:.TF1 (Log Double):.TL1 (Log Double):.EB (Log Double)
      TW (ITbl _ _ _ pf) _ = sZ
      bs' = assocs pf
      pssum = (Numeric.Log.sum $ Prelude.map snd bs') / (fromIntegral n - 1)
      bs = Prelude.map (second (/pssum)) bs'
  in bs
{-# NoInline edgeProbPartFun #-}

-- | Turn the edge probabilities into a score matrix.

edgeProbScoreMat :: [(EdgeBoundary I, Log Double)] -> ScoreMatrix (Log Double)
edgeProbScoreMat xs' = ScoreMatrix m V.empty V.empty
  where m = fromAssocs l h 0 xs
        l = (Z:.0:.0)
        h = (Z:.maximum [f | (f :-> _,_) <- xs']:.maximum [t | (_ :-> t,_) <- xs'])
        xs = [ ((Z:.f:.t),p) | (f :-> t, p) <- xs' ]

