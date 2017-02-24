
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

import qualified Data.Bijection.HashMap as B
import           ADP.Fusion.Core
import           ADP.Fusion.EdgeBoundary
import           ADP.Fusion.Set1
import           Data.PrimitiveArray hiding (toList)
import           Data.PrimitiveArray.ScoreMatrix
import           Diagrams.TwoD.ProbabilityGrid
import           FormalLanguage
import           ShortestPath.SHP.Grammar.EdgeProbIO
import           Data.Vector.Generic.Unstream

import           BioInf.MutationOrder.RNA
import           BioInf.MutationOrder.MinDist (ScaleFunction(..))



-- | Before using @aInside@ the @ScoreMat@ needs to be scaled
--
-- TODO the @Edge@ needs to be an @EdgeWithActive@ to get the active bits
-- on the left in the set.

aInside :: Monad m => ScaleFunction -> Landscape -> SigEdgeProb m (Log Double) (Log Double) (Int:.From:.To) (Int:.To)
aInside scaled Landscape{..} = SigEdgeProb
  { edge = \x (fset:.From f:.To t) ->
      let frna = rnas HM.! (BitSet fset)
          trna = rnas HM.! (BitSet fset `xor` bit t)
--          fene = centroidEnergy frna
--          tene = centroidEnergy trna
--          res' = scaleFunction scaled (tene - fene) / s
          res' = scaled frna trna
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Edge",(BitSet fset,f,t),frna,trna,' ',res',res,x,x*res) $
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
--          fene = centroidEnergy frna
--          tene = centroidEnergy trna
--          res' = scaleFunction scaled (tene - fene) / s
          res' = scaled frna trna
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Node",n,frna,trna,' ',res',res,x,res*x) $
#endif
          x * res
  , fini = \l (fset:.From f:.To t) r ->
      let frna = rnas HM.! (BitSet fset)
          trna = rnas HM.! (BitSet fset `xor` bit t)
--          fene = centroidEnergy frna
--          tene = centroidEnergy trna
--          res' = scaleFunction scaled (tene - fene) / s
          res' = scaled frna trna
          res  = Exp . negate $ res'
      in
#ifdef ADPFUSION_DEBUGOUTPUT
          traceShow ("Fini",(BitSet fset,f,t),frna,trna,' ',res',res,l,r,l*r*res) $
#endif
          l * r * res
  , h    = SM.foldl' (+) 0
--  , h    = \s -> do v :: V.Vector (Log Double) <- streamToVectorM s
--                    return $ Numeric.Log.sum v
  } where !n = fromIntegral mutationCount
{-# Inline aInside #-}



type TF1 x = TwITbl Id Unboxed EmptyOk (BS1 Last I)      x
type TL1 x = TwITbl Id Unboxed EmptyOk (BS1 Last O)      x
type EB  x = TwITbl Id Unboxed EmptyOk (EdgeBoundary C)   x



-- | Extract the individual partition scores.

edgeProbPartFun :: ScaleFunction -> Landscape -> ([(Boundary Last I, Log Double)], [(EdgeBoundary C, Log Double)])
edgeProbPartFun scaled landscape =
  let n       = mutationCount landscape
      (Z:.sF:.sL:.sZ) = mutateTablesST $ gEdgeProb (aInside scaled landscape)
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

edgeProbScoreMatrix :: Landscape -> [(EdgeBoundary C, Log Double)] -> ScoreMatrix (Log Double)
edgeProbScoreMatrix Landscape{..} xs' = ScoreMatrix m nprobs names names
  where m = fromAssocs l h 0 xs
        l = (Z:.0:.0)
        h = (Z:.maximum [f | (f :-> _,_) <- xs']:.maximum [t | (_ :-> t,_) <- xs'])
        (Z:._:.hh) = h
        xs = [ ((Z:.f:.t),p) | (f :-> t, p) <- xs' ]
        (_,Z:._:.n) = bounds m
        names = V.fromList [ T.pack . show . (+1)
                           . maybe (error "MutationOrder.EdgeProb.edgeProbScoreMatrix") id
                           $ B.lookupR mutationPositions k | k <- [0..n]
                           ]
        rowsums = HM.fromListWith (+) [ (f,p) | (f :-> t, p) <- xs' ]
        nprobs = fromAssocs 0 hh 1 [ (f,1-p) | f <- [0..hh], let p = rowsums HM.! f ]

