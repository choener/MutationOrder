
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
import           ShortestPath.SHP.MinDist

import           BioInf.MutationOrder.RNA



data ScaleFunction
  = ScaleId
  | ScalePositiveSquared
  deriving (Show,Data,Typeable)

-- | Minimal distance algebra
--
-- TODO The two Ints are the indices of the nodes and could be replaced?

aMinDist :: Monad m => ScaleFunction -> Landscape -> SigMinDist m Double Double (Int:.From:.To) Int
aMinDist scaled Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                      in  -- traceShow (BitSet fset, BitSet fset `setBit` f `setBit` t) $
                                          x + scaleFunction scaled (centroidEnergy trna - centroidEnergy frna)
  , mpty = \() -> 0
  , node = \n -> let frna = rnas HM.! (BitSet 0)
                     trna = rnas HM.! (BitSet 0 `setBit` n)
                 in  scaleFunction scaled $ centroidEnergy trna - centroidEnergy frna
  , fini = id
  , h    = SM.foldl' min 999999
  }
{-# Inline aMinDist #-}

scaleFunction ScaleId = id
scaleFunction ScalePositiveSquared = \x -> if x <= 0 then x else x*x
{-# Inline scaleFunction #-}

{-
-- | Maximum edge probability following the probabilities generated from
-- the @EdgeProb@ grammar.

aMaxEdgeProb :: Monad m => ScoreMat (Log Double) -> SigMinDist m (Log Double) (Log Double) (From:.To) Int
aMaxEdgeProb s = SigMinDist
  { edge = \x e -> x * (s .!. e)
  , mpty = \() -> 1
  , node = \n -> 1
  , fini = id
  , h    = SM.foldl' max 0
  }
{-# Inline aMaxEdgeProb #-}
-}

-- | This should give the correct order of nodes independent of the
-- underlying @Set1 First@ or @Set1 Last@ because the @(From:.To)@ system
-- is agnostic over these.
--
-- TODO Use text builder

aPretty :: Monad m => Landscape -> SigMinDist m Text [Text] (Int:.From:.To) Int
aPretty Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> let frna = rnas HM.! (BitSet fset)
                                          trna = rnas HM.! (BitSet fset `setBit` f `setBit` t)
                                          eM = mfeEnergy trna - mfeEnergy frna
                                          eC = centroidEnergy trna - centroidEnergy frna
                                          f' = fromJust $ B.lookupR mutationPositions f
                                          t' = fromJust $ B.lookupR mutationPositions t
                                      in  T.concat [x, showMut frna trna f' eM eC]
  , mpty = \()  -> ""
  , node = \n   -> let frna = rnas HM.! (BitSet 0)
                       trna = rnas HM.! (BitSet 0 `setBit` n)
                       n'   = fromJust $ B.lookupR mutationPositions n
                       eM   = mfeEnergy trna - mfeEnergy frna
                       eC   = centroidEnergy trna - centroidEnergy frna
                   in  T.concat [showHdr frna n', showMut frna trna n' eM eC]
  , fini = id
  , h    = SM.toList
  } where
      showHdr frna n = T.concat
        [ T.pack $ printf "mutation         mfe     centr   "
        , T.pack $ VU.toList $ VU.replicate (BS.length $ primarySequence frna) ' ' VU.// (map (,'v') . sort . map fst $ B.toList mutationPositions)
        , T.pack $ "\n" ++ replicate 33 ' '
        , T.pack . take (BS.length $ primarySequence frna) . concat $ zipWith (\xs x -> xs ++ show x) (repeat $ "    .    ") (drop 1 $ cycle [0..9])
        , "\n"
        , T.pack $ printf "ancestral        %5.1f   %5.1f   " (mfeEnergy frna) (centroidEnergy frna)
        , T.pack $ BS.unpack $ primarySequence frna
        , "\n"
        ]
      showMut frna trna n eM eC = T.concat
        [ T.pack $ printf "%5d            %5.1f   %5.1f   " (n+1) eM eC
        , T.pack . BS.unpack $ primarySequence trna
        , "\n"
        ]
{-# Inline aPretty #-}

-- | Count co-optimals

aCount :: Monad m => Landscape -> SigMinDist m Integer [Integer] (Int:.From:.To) Int
aCount Landscape{..} = SigMinDist
  { edge = \x (fset:.From f:.To t) -> x
  , mpty = \()  -> 1
  , node = \n   -> 1
  , fini = id
  , h    = \xs -> SM.foldl' (+) 0 xs >>= \x -> return [x]
  }
{-# Inline aCount #-}

{-
-- | Before using @aInside@ the @ScoreMat@ needs to be scaled
-- appropriately! Due to performance reasons we don't want to do this
-- within @aInside@.

aInside :: Monad m => ScoreMat (Log Double) -> SigMinDist m (Log Double) (Log Double) (From:.To) Int
aInside s = SigMinDist
  { edge = \x e -> s .!. e * x
  , mpty = \() -> 1
  , node = \n -> 1
  , fini = id
  , h    = SM.foldl' (+) 0
  }
{-# Inline aInside #-}
-}


type TS1 x = TwITbl Id Unboxed EmptyOk (BS1 First I)      x
type U   x = TwITbl Id Unboxed EmptyOk (Unit I)           x
type PF  x = TwITbl Id Unboxed EmptyOk (Boundary First I) x

type BT1 x b = TwITblBt Unboxed EmptyOk (BS1 First I) x Id Id b
type BTU x b = TwITblBt Unboxed EmptyOk (Unit I)      x Id Id b



-- | Run the minimal distance algebra.
--
-- This produces one-boundary sets. Meaning that for each boundary we get
-- the total distance within the set.

forwardMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1 Double:.U Double
forwardMinDist1 scaleFunction landscape =
  let n = mutationCount landscape
  in  mutateTablesST $ gMinDist (aMinDist scaleFunction landscape)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) (999999) []))
        (ITbl 1 0 EmptyOk (fromAssocs Unit         Unit                           (999999) []))
        EdgeWithSet
        Singleton
{-# NoInline forwardMinDist1 #-}

backtrackMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1 Double:.U Double -> [Text]
backtrackMinDist1 scaleFunction landscape (Z:.ts1:.u) = unId $ axiom b
  where !(Z:.bt1:.b) = gMinDist (aMinDist scaleFunction landscape <|| aPretty landscape)
                            (toBacktrack ts1 (undefined :: Id a -> Id a))
                            (toBacktrack u   (undefined :: Id a -> Id a))
                            EdgeWithSet
                            Singleton
                        :: Z:.BT1 Double Text:.BTU Double Text
{-# NoInline backtrackMinDist1 #-}

countBackMinDist1 :: ScaleFunction -> Landscape -> Z:.TS1 Double:.U Double -> [Integer]
countBackMinDist1 scaleFunction landscape (Z:.ts1:.u) = unId $ axiom b
  where !(Z:.bt1:.b) = gMinDist (aMinDist scaleFunction landscape <|| aCount landscape)
                            (toBacktrack ts1 (undefined :: Id a -> Id a))
                            (toBacktrack u   (undefined :: Id a -> Id a))
                            EdgeWithSet
                            Singleton
                        :: Z:.BT1 Double Integer:.BTU Double Integer
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

{-
-- | Extract the individual partition scores.

boundaryPartFun :: Double -> ScoreMat Double -> [(Boundary First I,Log Double)]
boundaryPartFun temperature landscape =
  let n       = numNodes landscape
      partMat = toPartMatrix temperature landscape
      (Z:.sM:.bM) = mutateTablesST $ gMinDist (aInside partMat)
                      (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) (-999999) []))
                      (ITbl 1 0 EmptyOk (fromAssocs (Boundary 0) (Boundary $ n-1)               (-999999) []))
                      Edge
                      Singleton
                    :: Z:.TS1 (Log Double):.PF (Log Double)
      TW (ITbl _ _ _ pf) _ = bM
      bs' = assocs pf
      pssum = Numeric.Log.sum $ Prelude.map snd bs'
      bs = Prelude.map (second (/pssum)) bs'
  in bs

{-# NoInline boundaryPartFun #-}

-- | Run the maximal edge probability grammar.

forwardMaxEdgeProb :: ScoreMat (Log Double) -> Z:.TS1 (Log Double):.U (Log Double)
forwardMaxEdgeProb landscape =
  let n = numNodes landscape
  in  mutateTablesST $ gMinDist (aMaxEdgeProb landscape)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))
        (ITbl 1 0 EmptyOk (fromAssocs Unit         Unit                           0 []))
        Edge
        Singleton
{-# NoInline forwardMaxEdgeProb #-}

backtrackMaxEdgeProb :: ScoreMat (Log Double) -> Z:.TS1 (Log Double):.U (Log Double) -> [Text]
backtrackMaxEdgeProb landscape (Z:.ts1:.u) = unId $ axiom b
  where !(Z:.bt1:.b) = gMinDist (aMaxEdgeProb landscape <|| aPretty landscape)
                            (toBacktrack ts1 (undefined :: Id a -> Id a))
                            (toBacktrack u   (undefined :: Id a -> Id a))
                            Edge
                            Singleton
                        :: Z:.BT1 (Log Double) Text:.BTU (Log Double) Text
{-# NoInline backtrackMaxEdgeProb #-}

-- | Given the @Set1@ produced in @forwardMinDist1@ we can now extract the
-- co-optimal paths using the @Set1 -> ()@ index change.
--
-- TODO do we want this one explicitly or make life easy and just extract
-- from all @forwardMinDist1@ paths?

runMaxEdgeProb :: ScoreMat (Log Double) -> (Log Double,[Text])
runMaxEdgeProb landscape = (unId $ axiom fwdu,bs)
  where !(Z:.fwd1:.fwdu) = forwardMaxEdgeProb landscape
        bs = backtrackMaxEdgeProb landscape (Z:.fwd1:.fwdu)
{-# NoInline runMaxEdgeProb #-}
-}

{-
test t fp = do
  sMat <- fromFile fp
  let (d,bt) = runCoOptDist sMat
  let ps = boundaryPartFun t sMat
  print d
  mapM_ print $ bt
  print $ length bt
  print $ length $ nub $ sort bt
  forM_ ps $ \(b,_) -> printf "%5s  " (sMat `nameOf` getBoundary b)
  putStrLn ""
  forM_ ps $ \(_,Exp p) -> printf "%0.3f  " (exp p)
  putStrLn ""
-}

