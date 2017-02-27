
-- | Calculate minimum-distance Hamiltonian Shortest Paths and
-- probabilities for starting nodes.
--
-- NOTE: We explicitly model starting nodes. For symmetrical distance
-- matrices, this reports begin/end probabilities. For asymmetrical
-- distance matrices, a second instances with @Last@ instead of @First@
-- boundary should be created to calculate begin/end probabilities
-- separately.

module BioInf.MutationOrder.MinDist where

import           Control.Arrow (second)
import           Control.Monad (forM_)
import           Data.Bits
import           Data.Data (Data)
import           Data.List (nub,sort)
import           Data.Text (Text)
import           Data.Typeable(Typeable)
import           Debug.Trace
import           Numeric.Log
import qualified Data.ByteString.Char8 as BS
import qualified Data.HashMap.Strict as HM
import qualified Data.Text as T
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import           Text.Printf
import qualified Data.Vector.Unboxed as VU
import           Data.Maybe (fromJust)

import qualified Data.Bijection.HashMap as B
import           ADP.Fusion.Core
import           ADP.Fusion.Set1
import           ADP.Fusion.Unit
import           Data.PrimitiveArray hiding (toList,map)
import           FormalLanguage
import           ShortestPath.SHP.Grammar.MinDist

import           BioInf.MutationOrder.RNA



-- | Given the 'RNA' we come from and the 'RNA' we mutate into, derive the
-- gain or loss by a scaling function.

type ScaleFunction = RNA -> RNA -> Double

-- | Minimal distance algebra
--
-- TODO The two Ints are the indices of the nodes and could be replaced?

aMinDist :: Monad m => ScaleFunction -> Landscape -> SigMinDist m Double Double (Int:.From:.To) (Int:.To)
aMinDist scaled Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                      in  -- traceShow (BitSet fset, BitSet fset `setBit` f `setBit` t) $
                                          -- x + scaleFunction scaled (centroidEnergy trna - centroidEnergy frna)
                                          x + scaled frna trna
  , mpty = \() -> 0
  , node = \(nset:.To n) ->
      let frna = rnas HM.! (BitSet 0)
          trna = rnas HM.! (BitSet 0 `setBit` n)
      in  -- scaleFunction scaled $ centroidEnergy trna - centroidEnergy frna
          scaled frna trna
  , fini = id
  , h    = SM.foldl' min 999999
  }
{-# Inline aMinDist #-}

-- | Sum over all states and collapse into boundary unscaled weights.

aInside :: Monad m => ScaleFunction -> Landscape -> SigMinDist m (Log Double) (Log Double) (Int:.From:.To) (Int:.To)
aInside scaled Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `xor` bit t)
                                      in
                                          -- traceShow ("edge",fset,f,t) $
                                          x + (Exp . negate $ scaled frna trna)
  , mpty = \() -> 1
  , node = \(nset:.To n) ->
      let frna = rnas HM.! (BitSet 0)
          trna = rnas HM.! (BitSet 0 `xor` bit n)
      in
          -- traceShow ("node",nset,n) .
          Exp . negate $ scaled frna trna
  , fini = id
  , h    = SM.foldl' (+) 0
  }
{-# Inline aInside #-}

-- | This should give the correct order of nodes independent of the
-- underlying @Set1 First@ or @Set1 Last@ because the @(From:.To)@ system
-- is agnostic over these.
--
-- TODO Use text builder

aPretty :: Monad m => ScaleFunction -> Landscape -> SigMinDist m Text [Text] (Int:.From:.To) (Int:.To)
aPretty scaled Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                          eM = mfeEnergy trna - mfeEnergy frna
                                          eC = centroidEnergy trna - centroidEnergy frna
                                          eS = scaled frna trna
                                          f' = fromJust $ B.lookupR mutationPositions f
                                          t' = fromJust $ B.lookupR mutationPositions t
                                      in  T.concat [x, showMut frna trna t' eM eC eS]
  , mpty = \()  -> ""
  , node = \(nset:.To n)  ->
      let
        frna = rnas HM.! (BitSet 0)
        trna = rnas HM.! (BitSet 0 `setBit` n)
        n'   = fromJust $ B.lookupR mutationPositions n
        eM   = mfeEnergy trna - mfeEnergy frna
        eC   = centroidEnergy trna - centroidEnergy frna
        eS   = scaled frna trna
      in  T.concat [showHdr frna n', showMut frna trna n' eM eC eS]
  , fini = id
  , h    = SM.toList
  } where
      showHdr frna n = T.concat
        [ T.pack $ printf "mutation         mfe    centr  scfun  "
        , T.pack $ VU.toList $ VU.replicate (BS.length $ primarySequence frna) ' ' VU.// (map (,'v') . sort . map fst $ B.toList mutationPositions)
        , T.pack $ "\n" ++ replicate 38 ' '
        , T.pack . take (BS.length $ primarySequence frna) . concat $ zipWith (\xs x -> xs ++ show x) (repeat $ "    .    ") (drop 1 $ cycle [0..9])
        , "\n"
        , T.pack $ printf "ancestral        %5.1f  %5.1f         " (mfeEnergy frna) (centroidEnergy frna)
        , T.pack $ BS.unpack $ primarySequence frna
        , "\n"
        ]
      showMut frna trna n eM eC eS = T.concat
        [ T.pack $ printf "%5d            %5.1f  %5.1f  %5.1f  " (n+1) eM eC eS
        , T.pack . BS.unpack $ primarySequence trna
        , "\n"
        ]
{-# Inline aPretty #-}

-- | Count co-optimals

aCount :: Monad m => Landscape -> SigMinDist m Integer [Integer] (Int:.From:.To) (Int:.To)
aCount Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> x
  , mpty = \()  -> 1
  , node = \n   -> 1
  , fini = id
  , h    = \xs -> SM.foldl' (+) 0 xs >>= \x -> return [x]
  }
{-# Inline aCount #-}



type TS1 x = TwITbl Id Unboxed EmptyOk (BS1 First I)      x
type U   x = TwITbl Id Unboxed EmptyOk (Unit I)           x
type PF  x = TwITbl Id Unboxed EmptyOk (Boundary First I) x

type TS1L x = TwITbl Id Unboxed EmptyOk (BS1 Last I)      x
type PFL  x = TwITbl Id Unboxed EmptyOk (Boundary Last I) x

type BT1 x b = TwITblBt Unboxed EmptyOk (BS1 First I) x Id Id b
type BTU x b = TwITblBt Unboxed EmptyOk (Unit I)      x Id Id b

type BT1L x b = TwITblBt Unboxed EmptyOk (BS1 Last I) x Id Id b



-- | Run the minimal distance algebra.
--
-- This produces one-boundary sets. Meaning that for each boundary we get
-- the total distance within the set.

forwardMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1L Double:.U Double
forwardMinDist1 scaleFunction landscape =
  let n = mutationCount landscape
  in  mutateTablesST $ gMinDist (aMinDist scaleFunction landscape)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) (999999) []))
        (ITbl 1 0 EmptyOk (fromAssocs Unit         Unit                           (999999) []))
        EdgeWithSet
        Singleton
{-# NoInline forwardMinDist1 #-}

backtrackMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1L Double:.U Double -> [Text]
backtrackMinDist1 scaleFunction landscape (Z:.ts1:.u) = unId $ axiom b
  where !(Z:.bt1:.b) = gMinDist (aMinDist scaleFunction landscape <|| aPretty scaleFunction landscape)
                            (toBacktrack ts1 (undefined :: Id a -> Id a))
                            (toBacktrack u   (undefined :: Id a -> Id a))
                            EdgeWithSet
                            Singleton
                        :: Z:.BT1L Double Text:.BTU Double Text
{-# NoInline backtrackMinDist1 #-}

countBackMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1L Double:.U Double -> [Integer]
countBackMinDist1 scaleFunction landscape (Z:.ts1:.u) = unId $ axiom b
  where !(Z:.bt1:.b) = gMinDist (aMinDist scaleFunction landscape <|| aCount landscape)
                            (toBacktrack ts1 (undefined :: Id a -> Id a))
                            (toBacktrack u   (undefined :: Id a -> Id a))
                            EdgeWithSet
                            Singleton
                        :: Z:.BT1L Double Integer:.BTU Double Integer
{-# NoInline countBackMinDist1 #-}

-- | Given the @Set1@ produced in @forwardMinDist1@ we can now extract the
-- co-optimal paths using the @Set1 -> ()@ index change.
--
-- TODO do we want this one explicitly or make life easy and just extract
-- from all @forwardMinDist1@ paths?

runCoOptDist :: ScaleFunction -> Landscape -> (Double,[Text])
runCoOptDist scaleFunction landscape = (unId $ axiom fwdu,bs)
  where !(Z:.fwd1:.fwdu) = forwardMinDist1 scaleFunction landscape
        bs = backtrackMinDist1 scaleFunction landscape (Z:.fwd1:.fwdu)
{-# NoInline runCoOptDist #-}

runCount :: ScaleFunction -> Landscape -> (Double,[Integer])
runCount scaleFunction landscape = (unId $ axiom fwdu,bs)
  where !(Z:.fwd1:.fwdu) = forwardMinDist1 scaleFunction landscape
        bs = countBackMinDist1 scaleFunction landscape (Z:.fwd1:.fwdu)
{-# NoInline runCount #-}

-- | Extract the individual partition scores.

boundaryPartFunFirst :: ScaleFunction -> Landscape -> [(Boundary First I,Log Double)]
boundaryPartFunFirst scaleFunction landscape =
  let n       = mutationCount landscape
      (Z:.sM:.bM) = mutateTablesST $ gMinDist (aInside scaleFunction landscape)
                      (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) (-999999) []))
                      (ITbl 1 0 EmptyOk (fromAssocs (Boundary 0) (Boundary $ n-1)               (-999999) []))
                      EdgeWithSet
                      Singleton
                    :: Z:.TS1 (Log Double):.PF (Log Double)
      TW (ITbl _ _ _ pf) _ = bM
      bs' = assocs pf
      pssum = Numeric.Log.sum $ Prelude.map snd bs'
      bs = Prelude.map (second (/pssum)) bs'
  in bs
{-# NoInline boundaryPartFunFirst #-}

boundaryPartFunLast :: ScaleFunction -> Landscape -> [(Boundary Last I,Log Double)]
boundaryPartFunLast scaleFunction landscape =
  let n       = mutationCount landscape
      (Z:.sM:.bM) = mutateTablesST $ gMinDist (aInside scaleFunction landscape)
                      (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) (-999999) []))
                      (ITbl 1 0 EmptyOk (fromAssocs (Boundary 0) (Boundary $ n-1)               (-999999) []))
                      EdgeWithSet
                      Singleton
                    :: Z:.TS1L (Log Double):.PFL (Log Double)
      TW (ITbl _ _ _ pf) _ = bM
      bs' = assocs pf
      pssum = Numeric.Log.sum $ Prelude.map snd bs'
      bs = Prelude.map (second (/pssum)) bs'
  in bs
{-# NoInline boundaryPartFunLast #-}

