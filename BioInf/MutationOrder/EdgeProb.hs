
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
import           ShortestPath.SHP.EdgeProbIO
import           Data.Vector.Generic.Unstream

import           BioInf.MutationOrder.RNA
import           BioInf.MutationOrder.MinDist (ScaleFunction(..),scaleFunction)



-- | Before using @aInside@ the @ScoreMat@ needs to be scaled
--
-- TODO the @Edge@ needs to be an @EdgeWithActive@ to get the active bits
-- on the left in the set.

aInside :: Monad m => ScaleFunction -> Landscape -> Double -> SigEdgeProb m (Log Double) (Log Double) (Int:.From:.To) (Int:.To)
aInside scaled Landscape{..} temperature = SigEdgeProb
  { edge = \x (fset:.From f:.To t) ->
      let frna = rnas HM.! (BitSet fset)
          trna = rnas HM.! (BitSet fset `xor` bit t)
          fene = centroidEnergy frna
          tene = centroidEnergy trna
          res' = scaleFunction scaled (tene - fene) / s
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Edge",(BitSet fset,f,t),frna,fene,trna,tene,' ',res',res,x,x*res) $
#endif
          x * res
  , mpty = \() ->
#ifdef ADPFUSION_DEBUGOUTPUT
                  traceShow "empty" $
#endif
                  1
  , node = \x (nset:.To n) ->
      let frna = rnas HM.! (BitSet nset)
          trna = rnas HM.! (BitSet nset `xor` bit n)
          fene = centroidEnergy frna
          tene = centroidEnergy trna
          res' = scaleFunction scaled (tene - fene) / s
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Node",n,frna,fene,trna,tene,' ',res',res,x,res*x) $
#endif
          x * res
  , fini = \l (fset:.From f:.To t) r ->
      let frna = rnas HM.! (BitSet fset)
          trna = rnas HM.! (BitSet fset `xor` bit t)
          fene = centroidEnergy frna
          tene = centroidEnergy trna
          res' = scaleFunction scaled (tene - fene) / s
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Fini",(BitSet fset,f,t),frna,fene,trna,tene,' ',res',res,l,r,l*r*res) $
#endif
          l * r * res
  , h    = SM.foldl' (+) 0
--  , h    = \s -> do v :: V.Vector (Log Double) <- streamToVectorM s
--                    return $ Numeric.Log.sum v
  } where !s = temperature
          !n = fromIntegral mutationCount
{-# Inline aInside #-}



type TF1 x = TwITbl Id Unboxed EmptyOk (BS1 Last I)      x
type TL1 x = TwITbl Id Unboxed EmptyOk (BS1 Last O)      x
type EB  x = TwITbl Id Unboxed EmptyOk (EdgeBoundary C)   x

{-
type BF1 x b = TwITblBt Unboxed EmptyOk (BS1 Last I)    x Id Id b
type BL1 x b = TwITblBt Unboxed EmptyOk (BS1 Last O)    x Id Id b
type BEB x b = TwITblBt Unboxed EmptyOk (EdgeBoundary I) x Id Id b
-}



-- | Extract the individual partition scores.

edgeProbPartFun :: ScaleFunction -> Double -> Landscape -> ([(Boundary Last I, Log Double)], [(EdgeBoundary C, Log Double)])
edgeProbPartFun scaled temperature landscape =
  let n       = mutationCount landscape
      (Z:.sF:.sL:.sZ) = mutateTablesST $ gEdgeProb (aInside scaled landscape temperature)
                          (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 1 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
                          (ITbl 2 0 EmptyOk (fromAssocs (0 :-> 0)    (0 :-> (n-1))                  0 []))
                          EdgeWithSet
                          Singleton
                        :: Z:.TF1 (Log Double):.TL1 (Log Double):.EB (Log Double)
      TW (ITbl _ _ _ pf ) _ = sZ
      TW (ITbl _ _ _ lkF) _ = sF
      TW (ITbl _ _ _ lkL) _ = sL
      bs' = assocs pf
      pssum  = (Numeric.Log.sum $ Prelude.map snd bs') / (fromIntegral n - 1)
      ibs'   = [ (Boundary b, p) | b <- [0..n] , let p = lkF ! (BS1 (2^n-1) (Boundary b)) ]
      obs'   = [ (Boundary b, p) | b <- [0..n] , let p = lkF ! (BS1 (2^n-1) (Boundary b)) ]
      ibssum = Numeric.Log.sum $ Prelude.map snd ibs'
      obssum = Numeric.Log.sum $ Prelude.map snd ibs'
      ibs    = Prelude.map (second (/ibssum)) ibs'
      bs     = Prelude.map (second (/ibssum)) bs'
  in
#ifdef ADPFUSION_DEBUGOUTPUT
      traceShow (assocs lkF) $
      traceShow (assocs lkL) $
      traceShow (assocs pf) $
      traceShow (ibssum,obssum) $
      traceShow (bs',pssum,bs) $
#endif
      (ibs,bs)
{-# NoInline edgeProbPartFun #-}

-- | Turn the edge probabilities into a score matrix.

edgeProbScoreMat :: [(EdgeBoundary I, Log Double)] -> ScoreMatrix (Log Double)
edgeProbScoreMat xs' = ScoreMatrix m V.empty V.empty
  where m = fromAssocs l h 0 xs
        l = (Z:.0:.0)
        h = (Z:.maximum [f | (f :-> _,_) <- xs']:.maximum [t | (_ :-> t,_) <- xs'])
        xs = [ ((Z:.f:.t),p) | (f :-> t, p) <- xs' ]

