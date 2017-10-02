
-- | Specialized mutation order grammar and algebras that incorporate exactly
-- one intermediate mutation. Not necessarily only of the backmutation kind.
--
-- Calculated are (i) the total evidence @Z@ for two input sequences. (ii) The
-- minimal weight distance.

module BioInf.MutationOrder.BackMutations where

import           Data.Bits
import           Data.HashMap.Strict (HashMap)
import qualified Data.HashMap.Strict as HM
import qualified Data.Vector.Fusion.Stream.Monadic as SM

import           Data.PrimitiveArray hiding (toList,map)
import           FormalLanguage
import           ADP.Fusion.Core hiding (PeekIndex)
import           ADP.Fusion.Set1 hiding (PeekIndex)
import           ADP.Fusion.Term.PeekIndex.Set1
import           ADP.Fusion.Unit hiding (PeekIndex)

import           BioInf.MutationOrder.SequenceDB (RNAfoldResult)

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



type ScaleFunction = RNAfoldResult → RNAfoldResult → Double

-- | Minimal distance calculation under the influence of one intermediate
-- mutation.
--
-- NOTE be very careful to check that @rnas@ (no intermediate mutation) and
-- @ntrs@ are used correctly!

aMinDist
  ∷ Monad m
  ⇒ Double
  -- ^ neutral element for @omin@.
  → (Double → Double → Double)
  -- ^ omin
  → (Double → Double → Double)
  -- ^ oplus
  → ScaleFunction
  -- ^ Combine two 'RNAfoldResult's into the resulting @Double@ score.
  → Int
  -- ^ Index position of the intermediate mutation or @(-1)@ if independent of the observed mutations.
  → HashMap Int RNAfoldResult
  -- ^ RNAs without the intermediate mutation set
  → HashMap Int RNAfoldResult
  -- ^ RNAs with the intermediate mutation set
  → SigMinDist m Double Double (Int:.From:.To) (Int:.To) (BS1 First I)
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
                                              trna = rnas HM.! (fset `setBit` f `setBit` t)
                                          in  if ipos < 0 || ipos /= t
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
                                              tset = fset `setBit` t
                                          in  if ipos < 0 || (not $ tset `testBit` ipos)
                                                then x `oplus` scaled frna trna
                                                else neutral
  -- ^ These are edges inserted while the intermediate mutation is active. Only
  -- allow intermediate mutations that do not move the @ipos@ mutation into its
  -- final state.
  -- TODO check if @f@ or @t@ will be set here
  , finalFromInter = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                                   trna = rnas HM.! tset
                                                   tset = fset `setBit` t
                                               in  if ipos >= 0 && ipos == t
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
                                              tset = fset `setBit` t
                                          in  if ipos < 0 || (fset `testBit` ipos)
                                                then x `oplus` scaled frna trna
                                                else neutral
  , finis = id
  -- ^ Collapse over possible end points
  , h = SM.foldl' omin neutral
  -- ^ Find the minimal distance from ancestral to extant sequence.
  }
{-# Inline aMinDist #-}

type FwdBS1  x = TwITbl Id Unboxed EmptyOk (BS1 First I) x
type FwdUnit x = TwITbl Id Unboxed EmptyOk (Unit      I) x

-- |
--
-- TODO check if inlining of @ScaleFunction@ improves performance
-- substantially.

forwardMinDist
  ∷ Int
  → ScaleFunction
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

-- |
--
-- TODO check if inlining of @ScaleFunction@ improves performance
-- substantially.

forwardEvidence
  ∷ Int
  → ScaleFunction
  → Int
  → HashMap Int RNAfoldResult
  → HashMap Int RNAfoldResult
  → Z:.FwdBS1 Double:.FwdBS1 Double:.FwdBS1 Double:.FwdUnit Double
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



