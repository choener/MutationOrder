
-- | Specialized mutation order grammar and algebras that incorporate exactly
-- one intermediate mutation. Not necessarily only of the backmutation kind.
--
-- Calculated are (i) the total evidence @Z@ for two input sequences. (ii) The
-- minimal weight distance.

module BioInf.MutationOrder.BackMutations where

import           ADP.Fusion.Core
import           Data.Bits
import           Data.HashMap.Strict (HashMap)
import           Data.PrimitiveArray hiding (toList,map)
import           FormalLanguage
import qualified Data.HashMap.Strict as HM
import qualified Data.Vector.Fusion.Stream.Monadic as SM

import           BioInf.MutationOrder.SequenceDB (RNAfoldResult)

[formalLanguage|
Verbose
Grammar: MinDist
N: X
N: B
N: U
N: S
T: n -- single node
T: k -- normal edge
T: b -- insert backmutation
T: u -- undo backmutation
S: S
-- A number of mutations that do not include the intermediate mutation.
X -> mpty <<< ε
X -> node <<< n
X -> edge <<< X k
-- Insert the single intermediate mutation, followed by more known events.
B -> edgB <<< X b
B -> edgI <<< B k
-- Undo the intermediate mutation, followed by more known events.
U -> edgU <<< B u
U -> edge <<< U k
S -> fini <<< U
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
  ⇒ ScaleFunction
  -- ^ Combine two 'RNAfoldResult's into the resulting @Double@ score.
  → HashMap Int RNAfoldResult
  -- ^ RNAs without the intermediate mutation set
  → HashMap Int RNAfoldResult
  -- ^ RNAs with the intermediate mutation set
  → SigMinDist m Double Double (Int:.From:.To) (Int:.From:.To) (Int:.To) (Int:.From:.To)
aMinDist scaled rnas ntrs = SigMinDist
  { edgB = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                         trna = ntrs HM.! fset
                                     in  x + scaled frna trna
  -- ^ The intermediate mutation is about to be set.
  , edgI = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                         trna = ntrs HM.! (fset `setBit` f `setBit` t)
                                     in  x + scaled frna trna
  -- ^ These are edges inserted while the intermediate mutation is active.
  , edgU = \x (fset:.From f:.To t) → let frna = ntrs HM.! fset
                                         trna = rnas HM.! fset
                                     in  x + scaled frna trna
  -- ^ Now the intermediate mutation is undone.
  , edge = \x (fset:.From f:.To t) → let frna = rnas HM.! fset
                                         trna = rnas HM.! (fset `setBit` f `setBit` t)
                                     in  x + scaled frna trna
  -- ^ These edges are mutational events without the influence of the
  -- intermediate mutation present.
  , fini = id
  -- ^ Collapse over possible end points
  , mpty = \() → 0
  -- ^ Not a single mutation has happened.
  , node = \(nset:.To n) → let frna = rnas HM.! 0
                               trna = rnas HM.! (0 `setBit` n)
                           in  scaled frna trna
  -- ^ Activate a single, first mutation.
  , h    = SM.foldl' min 999999
  -- ^ Find the minimal distance from ancestral to extant sequence.
  }
{-# Inline aMinDist #-}

type FwdBS1  x = TwITbl Id Unboxed EmptyOk (BS1 First I) x
type FwdUnit x = TwITbl Id Unboxed EmptyOk (Unit      I) x

forwardMinDist
  ∷ Int
  → ScaleFunction
  → HashMap Int RNAfoldResult
  → HashMap Int RNAfoldResult
  → Z:.FwdBS1 Double:.FwdUnit Double:.FwdBS1 Double:.FwdBS1 Double
forwardMinDist n scaled rnas ntrs =
  let 
  in  mutateTablesST $ gMinDist (aMinDist scaled rnas ntrs)
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))
        (ITbl 0 0 EmptyOk (fromAssocs Unit         Unit                           999999 []))
        (ITbl 0 0 EmptyOk (fromAssocs (BS1 0 (-1)) (BS1 (2^n-1) (Boundary $ n-1)) 999999 []))
{-# Inline forwardMinDist #-}

