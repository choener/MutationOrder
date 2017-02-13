
-- | Here we collect the necessary data structures for the RNAs to be
-- compared. This data is serialized to disk once calculated, since it is
-- most likely the part that takes longest.
--
-- TODO if the vienna wrapper allows, we should parallelize the
-- calculations.
--
-- TODO nice interface counting up?

module BioInf.MutationOrder.RNA where



-- | A single RNA with pre-calculated elements.
--
-- TODO include the basepair probability matrix? Can we "compress" that
-- one?

data RNA = RNA
  { mutationSet :: VU.Vector (Int,Char)
    -- ^ we store just the mutation set, since this is more sparse and
    -- gives access to the mutational events.
  , primarySequence :: ByteString
    -- ^ store RNA sequence too, for now
  , structure :: ByteString
    -- ^ the structure we get
  }
  deriving (Show,Eq,Generics)

data Landscape = Landscape
  { rnas :: V.Vector RNA
    -- ^ the individual RNA mutations. The index should be calculated from
    -- @linearIndex 0 high mutationSet@
  , mutationCount :: Int
    -- ^ how many nucleotides are mutated in total
  }
  deriving (Show,Eq,Generics)

createRNAlandscape :: ByteString -> ByteString -> Landscape
createRNAlandscape = origin mutation = Landscape
  { rnas = undefined
  }
  where

