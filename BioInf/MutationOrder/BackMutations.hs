
-- | Specialized mutation order grammar and algebras that incorporate exactly
-- one intermediate mutation. Not necessarily only of the backmutation kind.
--
-- Calculated are (i) the total evidence @Z@ for two input sequences. (ii) The
-- minimal weight distance.

module BioInf.MutationOrder.BackMutations where

import qualified Data.Char as C
import           Data.Bits
import           Data.HashMap.Strict (HashMap)
import qualified Data.HashMap.Strict as HM
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import           Numeric.Log
import           Data.Text (Text)
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as BS
import           Text.Printf
import qualified Data.Vector.Unboxed as VU

import qualified Data.Bijection.HashMap as B
import           Data.PrimitiveArray hiding (toList,map)
import qualified Data.PrimitiveArray as PA
import           FormalLanguage
import           ADP.Fusion.Core hiding (PeekIndex)
import           ADP.Fusion.Set1 hiding (PeekIndex)
import           ADP.Fusion.Term.PeekIndex.Set1
import           ADP.Fusion.Unit hiding (PeekIndex)

import           BioInf.MutationOrder.SequenceDB (RNAfoldResult(..))

[formalLanguage|
Verbose
Grammar: MinDist
N: Begin
N: Inter
N: Final
N: Start
T: peek -- insert backmutation
T: edg  -- normal edge
T: nd   -- single node
S: Start
-- A number of mutations that do not include the intermediate mutation.
Begin -> beginEmpty <<< ε
Begin -> beginNode  <<< nd
Begin -> beginEdge  <<< Begin edg
-- Insert the single intermediate mutation, followed by more known events.
-- TODO this *could* be done with a statically active guard!
Inter -> interFromBegin <<< Begin peek  -- activate backmutation or intermediate mutation
Inter -> interEdge      <<< Inter edg
-- Undo the intermediate mutation, followed by more known events.
Final -> finalFromInter   <<< Inter edg   -- this moves from an intermediate to a final mutation
Final -> finalUnBackmut   <<< Inter peek  -- deactivates the backmutation
Final -> finalEdge        <<< Final edg
Start -> finis            <<< Final
//
Emit: MinDist
|]

makeAlgebraProduct ''SigMinDist



type ScaleFunction d = RNAfoldResult → RNAfoldResult → d

-- | Minimal distance calculation under the influence of one intermediate
-- mutation.
--
-- NOTE be very careful to check that @rnas@ (no intermediate mutation) and
-- @ntrs@ are used correctly!

aMinDist
  ∷ Monad m
  ⇒ d
  -- ^ neutral element for @omin@.
  → (d → d → d)
  -- ^ omin
  → (d → d → d)
  -- ^ oplus
  → ScaleFunction d
  -- ^ Combine two 'RNAfoldResult's into the resulting @d@ score.
  → Int
  -- ^ Index position of the intermediate mutation or @(-1)@ if independent of the observed mutations.
  -- TODO IN COORDINATES OF THE BITSET SPACE!
  → HashMap Int RNAfoldResult
  -- ^ RNAs without the intermediate mutation set
  → HashMap Int RNAfoldResult
  -- ^ RNAs with the intermediate mutation set
  → SigMinDist m d d (Int:.From:.To) (Int:.To) (BS1 First I)
aMinDist neutral omin oplus scaled ipos rnas ntrs = SigMinDist
  { beginEmpty = \() → neutral
  -- ^ Not a single mutation has happened.
  , beginNode = \(nset:.To n) → let frna = rnas HM.! 0
                                    trna = rnas HM.! (0 `setBit` n)
                                in  if ipos < 0 || n /= ipos
                                      then scaled frna trna
                                      else neutral
  -- ^ Activate a single, first mutation.
  , beginEdge = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                              trna = rnas HM.! tset
                                              tset = fset `setBit` f `setBit` t
                                          in  if ipos < 0 || (not $ tset `testBit` ipos)
                                                then x `oplus` scaled frna trna
                                                else neutral
  -- ^ These edges are mutational events without the influence of the
  -- intermediate mutation present.
  , interFromBegin = \x (BS1 fset b) → let frna = rnas HM.! getBitSet fset
                                           trna = ntrs HM.! getBitSet fset
                                       -- Only allow if either the mutation does not
                                       -- influence the observed mutations or the
                                       -- influenced mutation is still available.
                                       in  if ipos < 0 || (not $ fset `testBit` ipos)
                                            then x `oplus` scaled frna trna
                                            else neutral
  -- ^ The intermediate mutation is about to be set. Only if @ipos@ is
  -- independent of the other mutational events or we have not tried activating it.
  , interEdge = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                              trna = ntrs HM.! tset
                                              tset = fset `setBit` f `setBit` t
                                          in  if ipos < 0 || (not $ tset `testBit` ipos)
                                                then x `oplus` scaled frna trna
                                                else neutral
  -- ^ These are edges inserted while the intermediate mutation is active. Only
  -- allow intermediate mutations that do not move the @ipos@ mutation into its
  -- final state.
  -- TODO check if @f@ or @t@ will be set here
  , finalFromInter = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                                   trna = rnas HM.! tset
                                                   tset = fset `setBit` f `setBit` t
                                               in  if ipos == f
                                                     then x `oplus` scaled frna trna
                                                     else neutral
  -- ^ Flips the intermediate mutation (@ipos >= 0@) into the final state.
  -- TODO check if @t@ is the one being set
  , finalUnBackmut = \x (BS1 fset t) → let frna = ntrs HM.! getBitSet fset
                                           trna = rnas HM.! getBitSet fset
                                       -- Only allow if either the mutation does not
                                       -- influence the observed mutations or the
                                       -- influenced mutation is 
                                       in  if ipos < 0
                                            then x `oplus` scaled frna trna
                                            else neutral
  -- ^ Now the intermediate mutation is undone. Is only used if @ipos@ is
  -- independent (@==(-1)@) of the observed events.
  , finalEdge = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                              trna = rnas HM.! tset
                                              tset = fset `setBit` f `setBit` t
                                          in  if ipos < 0 || (fset `testBit` ipos)
                                                then x `oplus` scaled frna trna
                                                else neutral
  , finis = id
  -- ^ Collapse over possible end points
  , h = SM.foldl' omin neutral
  -- ^ Find the minimal distance from ancestral to extant sequence.
  }
{-# Inline aMinDist #-}

-- |
--
-- TODO Use text builder

aPretty
  ∷ Monad m
  ⇒ ScaleFunction Double
  → Int
  -- ^ index of intermediate mutation
  → Int
  -- ^ real ipos value
  → HashMap Int RNAfoldResult
  -- ^ rnas
  → HashMap Int RNAfoldResult
  -- ^ intermediate rnas
  → B.BimapHashMap Int Int
  -- ^ actual mutation position / bit in bitset
  → SigMinDist m Text [Text] (Int:.From:.To) (Int:.To) (BS1 First I)
aPretty scaled ipos realPos rnas ntrs mutpos = SigMinDist
  { beginEmpty = \() → ""
  , beginNode = \(nset:.To n) → let frna = rnas HM.! 0
                                    trna = rnas HM.! (0 `setBit` n)
                                in  if ipos < 0 || n /= ipos
                                      then T.concat [showHdr frna, showMut frna trna n]
                                      else "BEGINNODE ERROR\n"
  -- ^ Activate a single, first mutation.
  , beginEdge = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                              trna = rnas HM.! (fset `setBit` f `setBit` t)
                                          in  if ipos < 0 || ipos /= f
                                                then T.concat [x, showMut frna trna f]
                                                else "BEGINEDGE ERROR\n"
  -- ^ These edges are mutational events without the influence of the
  -- intermediate mutation present.
  , interFromBegin = \x (BS1 fset b) → let frna = rnas HM.! getBitSet fset
                                           trna = ntrs HM.! getBitSet fset
                                       -- Only allow if either the mutation does not
                                       -- influence the observed mutations or the
                                       -- influenced mutation is still available.
                                       in  if ipos < 0 || (not $ fset `testBit` ipos)
                                            then T.concat [x, showBackmut frna trna realPos]
                                            else "INTERFROMBEGIN ERROR\n"
  -- ^ The intermediate mutation is about to be set. Only if @ipos@ is
  -- independent of the other mutational events or we have not tried activating it.
  , interEdge = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                              trna = ntrs HM.! tset
                                              tset = fset `setBit` f `setBit` t
                                          in  if ipos < 0 || (not $ tset `testBit` ipos)
                                                then T.concat [x, showMut frna trna f]
                                                else "INTEREDGE ERROR\n"
  -- ^ These are edges inserted while the intermediate mutation is active. Only
  -- allow intermediate mutations that do not move the @ipos@ mutation into its
  -- final state.
  -- TODO check if @f@ or @t@ will be set here
  , finalFromInter = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                                   trna = rnas HM.! tset
                                                   tset = fset `setBit` f `setBit` t
                                               in  if ipos == f
                                                     then T.concat [x, showUnBackmut frna trna realPos] -- , T.pack (show (ipos,BitSet fset,BitSet tset,f,t)), "\n"]
                                                     else "FINALFROMINTER ERROR\n"
  -- ^ Flips the intermediate mutation (@ipos >= 0@) into the final state.
  -- TODO check if @t@ is the one being set
  , finalUnBackmut = \x (BS1 fset t) → let frna = ntrs HM.! getBitSet fset
                                           trna = rnas HM.! getBitSet fset
                                       -- Only allow if either the mutation does not
                                       -- influence the observed mutations or the
                                       -- influenced mutation is 
                                       in  if ipos < 0
                                            then T.concat [x, showUnBackmut frna trna realPos]
                                            else "FINALUNBACKMUT ERROR\n"
  -- ^ Now the intermediate mutation is undone. Is only used if @ipos@ is
  -- independent (@==(-1)@) of the observed events.
  , finalEdge = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                              trna = rnas HM.! tset
                                              tset = fset `setBit` f `setBit` t
                                          in  if ipos < 0 || (fset `testBit` ipos)
                                                then T.concat [x, showMut frna trna f]
                                                else "FINALEDGE ERROR\n"
  , finis = id
  -- ^ Collapse over possible end points
  , h = SM.toList
  -- ^ Find the minimal distance from ancestral to extant sequence.
  } where
      muts = let as = VU.generate n (const ' ')
                 n  = BS.length $ rnaFoldSequence $ rnas HM.! 0
                 bs = as VU.// [ (k*10, C.chr $ C.ord '0' + k `mod` 10) | k ← [1 .. n `div` 10] ]
                 cs = bs VU.// [ (k, 'v') | (k,_) ← B.toList mutpos ]
                 ds = cs VU.// [ if ipos < 0 then (realPos-1,'!') else (realPos-1,'+') ]
             in  VU.toList ds -- ++ show (ipos,realPos)
      showHdr ∷ RNAfoldResult → Text
      showHdr frna = T.pack $ printf "%s\n            %s\n%5.1f       %s\n%5.1f       %s\n"
                                (replicate 12 ' ' ++ muts)
                                (BS.unpack $ rnaFoldSequence frna)
                                (rnaFoldMFEEner frna)
                                (BS.unpack $ rnaFoldMFEStruc frna)
                                (rnaFoldCentroidEner frna)
                                (BS.unpack $ rnaFoldCentroidStruc frna)
      showMut ∷ RNAfoldResult → RNAfoldResult → Int → Text
      showMut frna trna p = T.pack $ printf "%5d %5.1f %s\n%s%s\n%s%s\n"
                                        (maybe (error $ "showMut" ++ show p) (+1) $ B.lookupR mutpos p)
                                        (deltaE frna trna)
                                        (BS.unpack $ rnaFoldSequence trna)
                                        (replicate 12 ' ')
                                        (BS.unpack $ rnaFoldMFEStruc trna)
                                        (replicate 12 ' ')
                                        (BS.unpack $ rnaFoldCentroidStruc trna)
      showBackmut ∷ RNAfoldResult → RNAfoldResult → Int → Text
      showBackmut frna trna p = T.pack $ printf "%c %3d %5.1f %s\n%s%s\n%s%s\n"
                                        (if ipos < 0 then '!' else '+')
                                        p
                                        (deltaE frna trna)
                                        (BS.unpack $ rnaFoldSequence trna)
                                        (replicate 12 ' ')
                                        (BS.unpack $ rnaFoldMFEStruc trna)
                                        (replicate 12 ' ')
                                        (BS.unpack $ rnaFoldCentroidStruc trna)
      showUnBackmut ∷ RNAfoldResult → RNAfoldResult → Int → Text
      showUnBackmut frna trna p = T.pack $ printf "%c %3d %5.1f %s\n%s%s\n%s%s\n"
                                             (if ipos < 0 then '!' else '+')
                                             p
                                             (deltaE frna trna)
                                             (BS.unpack $ rnaFoldSequence trna)
                                             (replicate 12 ' ')
                                             (BS.unpack $ rnaFoldMFEStruc trna)
                                             (replicate 12 ' ')
                                             (BS.unpack $ rnaFoldCentroidStruc trna)
      deltaE ∷ ScaleFunction Double
      deltaE = scaled -- frna trna = rnaFoldMFEEner trna - rnaFoldMFEEner frna
{-# Inline aPretty #-}

type FwdBS1  x = TwITbl Id Unboxed EmptyOk (BS1 First I) x
type FwdUnit x = TwITbl Id Unboxed EmptyOk (Unit      I) x

type BTS1   x b = TwITblBt Unboxed EmptyOk (BS1 First I) x Id Id b
type BTUnit x b = TwITblBt Unboxed EmptyOk (Unit I)      x Id Id b

-- |
--
-- TODO check if inlining of @ScaleFunction@ improves performance
-- substantially.

forwardMinDist
  ∷ Int
  → ScaleFunction Double
  → Int
  → HashMap Int RNAfoldResult
  → HashMap Int RNAfoldResult
  → Z:.FwdBS1 Double:.FwdBS1 Double:.FwdBS1 Double:.FwdUnit Double
forwardMinDist n scaled ipos rnas ntrs =
  let 
  in  mutateTablesST $ gMinDist (aMinDist 999999 min (+) scaled ipos rnas ntrs)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))   -- Begin
        (ITbl 2 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))   -- Final
        (ITbl 1 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))   -- Inter
        (ITbl 3 0 EmptyOk (fromAssocs Unit         Unit                           999999 []))   -- Start
        EdgeWithSet   -- normal mutational event
        Singleton     -- first mutational event
        (PeekIndex ∷ PeekIndex (BS1 First I))     -- undo intermediate
{-# NoInline forwardMinDist #-}

forwardMinDistValue fwd = md PA.! PA.Unit
  where (Z:.fwdB:.fwdF:.fwdI:.fwdS) = fwd
        TW (ITbl _ _ _ md) _ = fwdS

backtrackMinDist1
  ∷ Int
  → ScaleFunction Double
  → Int
  → Int
  → HashMap Int RNAfoldResult
  → HashMap Int RNAfoldResult
  → B.BimapHashMap Int Int
  → Z:.FwdBS1 Double:.FwdBS1 Double:.FwdBS1 Double:.FwdUnit Double
  → [Text]
backtrackMinDist1 n scaled ipos realPos rna ntrs mutpos (Z:.fwdB:.fwdF:.fwdI:.fwdS) = unId $ axiom btS
  where !(Z:.btB:.btF:.btI:.btS) =
            gMinDist (aMinDist 999999 min (+) scaled ipos rna ntrs <|| aPretty scaled ipos realPos rna ntrs mutpos)
                            (toBacktrack fwdB (undefined :: Id a -> Id a))
                            (toBacktrack fwdF (undefined :: Id a -> Id a))
                            (toBacktrack fwdI (undefined :: Id a -> Id a))
                            (toBacktrack fwdS (undefined :: Id a -> Id a))
                            EdgeWithSet
                            Singleton
                            (PeekIndex ∷ PeekIndex (BS1 First I))
                        :: Z:.BTS1 Double Text:.BTS1 Double Text:.BTS1 Double Text:.BTUnit Double Text
{-# NoInline backtrackMinDist1 #-}

-- |
--
-- TODO check if inlining of @ScaleFunction@ improves performance
-- substantially.

forwardEvidence
  ∷ Int
  → ScaleFunction (Log Double)
  → Int
  → HashMap Int RNAfoldResult
  → HashMap Int RNAfoldResult
  → Z:.FwdBS1 (Log Double):.FwdBS1 (Log Double):.FwdBS1 (Log Double):.FwdUnit (Log Double)
forwardEvidence n scaled ipos rnas ntrs =
  let 
  in  mutateTablesST $ gMinDist (aMinDist 0 (+) (*) scaled ipos rnas ntrs)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))   -- Begin
        (ITbl 2 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))   -- Final
        (ITbl 1 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 0 []))   -- Inter
        (ITbl 3 0 EmptyOk (fromAssocs Unit         Unit                           0 []))   -- Start
        EdgeWithSet   -- normal mutational event
        Singleton     -- first mutational event
        (PeekIndex ∷ PeekIndex (BS1 First I))     -- undo intermediate
{-# NoInline forwardEvidence #-}

forwardEvidenceValue fwd = md PA.! PA.Unit
  where (Z:.fwdB:.fwdF:.fwdI:.fwdS) = fwd
        TW (ITbl _ _ _ md) _ = fwdS



