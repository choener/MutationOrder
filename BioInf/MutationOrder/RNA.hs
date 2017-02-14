
-- | Here we collect the necessary data structures for the RNAs to be
-- compared. This data is serialized to disk once calculated, since it is
-- most likely the part that takes longest.
--
-- TODO if the vienna wrapper allows, we should parallelize the
-- calculations.
--
-- TODO nice interface counting up?

module BioInf.MutationOrder.RNA where

import           Data.Bits
import           Codec.Compression.GZip (compress,decompress)
import           Control.Arrow (second)
import           Control.DeepSeq
import           Control.Parallel.Strategies
import           Data.ByteString (ByteString)
import           Data.Maybe (catMaybes)
import           Data.Serialize
import           Data.Vector.Serialize
import           Data.Vector.Strategies
import           Debug.Trace
import           GHC.Generics
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy as BSL
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import           System.IO.Unsafe (unsafePerformIO)
import qualified Data.HashMap.Strict as HM
import           Data.Serialize.Instances

import qualified Data.Bijection.HashMap as B
import           BioInf.ViennaRNA.Bindings
import qualified Data.PrimitiveArray as PA



-- | A single RNA with pre-calculated elements.
--
-- All calculations are done at 37 C.
--
-- TODO include the basepair probability matrix? Can we "compress" that
-- one?

data RNA = RNA
  { mutationSet       :: !(VU.Vector (Int,Char))
    -- ^ we store just the mutation set, since this is more sparse and
    -- gives access to the mutational events.
  , primarySequence   :: !ByteString
    -- ^ store RNA sequence too, for now
  , mfeStructure      :: !ByteString
    -- ^ the mfe structure we get
  , mfeEnergy         :: !Double
    -- ^ mfe energy of the structure
  , centroidStructure :: !ByteString
    -- ^ the centroid structure
  , centroidEnergy    :: !Double
  }
  deriving (Show,Eq,Generic)

instance NFData     RNA
instance Serialize  RNA

-- | Given the primary sequence and the mutation set, fill the 'RNA'
-- structure.
--
-- NOTE This wraps some @ViennaRNA-bindings@ calls that are in @IO@.
--
-- TODO check if these calls are *really* thread-safe!

mkRNA
  :: ByteString
  -- ^ primary sequence of the *origin* RNA
  -> VU.Vector (Int,Char)
  -- ^ set of mutations compared to the origin
  -> RNA
mkRNA inp' ms = RNA
  { mutationSet       = ms
  , primarySequence   = inp
  , mfeStructure      = s
  , mfeEnergy         = e
  , centroidStructure = BS.empty
  , centroidEnergy    = 0
  }
  where
    inp   = insertMutations ms inp'
    (e,s) = second BS.pack . unsafePerformIO . mfeTemp 37 $ BS.unpack inp
    (cE,cS) = second BS.pack . unsafePerformIO . centroidTemp 37 $ BS.unpack inp

-- | Insert a set of mutations in a @ByteString@.

insertMutations :: VU.Vector (Int,Char) -> ByteString -> ByteString
insertMutations ms s' = VU.foldl' go s' ms
  where go s (k,c) =
          let (h,t) = BS.splitAt k s
          in  BS.concat [h, BS.singleton c, BS.drop 1 t]

data Landscape = Landscape
  { rnas                  :: HM.HashMap (PA.BitSet PA.I) RNA
    -- ^ the individual RNA mutations. The index should be calculated from
    -- @linearIndex 0 high mutationSet@
  , mutationCount         :: !Int
    -- ^ how many nucleotides are mutated in total
  , landscapeOrigin       :: !ByteString
    -- ^ the ancestral sequence
  , landscapeDestination  :: !ByteString
    -- ^ the final sequence
  }
  deriving (Show,Eq,Generic)

instance NFData     Landscape
instance Serialize  Landscape

-- |
--
-- TODO prime candidate for parallelization. ViennaRNA-bindings currently
-- does not allow parallel runs! It would be possible to consider
-- externalizing this, but for now we just run single-threaded.

createRNAlandscape :: Bool -> ByteString -> ByteString -> Landscape
createRNAlandscape verbose origin mutation = Landscape
  { rnas                  = rs -- `using` (parVector chunkSize)
  , mutationCount         = length . filter (>1) . map length $ pms
  , landscapeOrigin       = origin
  , landscapeDestination  = mutation
  }
  where
    rs  = HM.fromList . map pairWithBitSet $ zipWith talk mus [0..]
    talk s c = (if (c `mod` 1000 == 0 && verbose) then traceShow c else id) mkRNA origin s
    mus = map (VU.fromList . catMaybes)
        . sequence
        $ pms
    -- possible mutations
    pms = zipWith3 genM (BS.unpack origin) (BS.unpack mutation) [0..]
    genM a b k | a==b = [Nothing]
               | otherwise = [Nothing,Just (k,b)]
    -- pair each @RNA@ with the correct bitset
    pairWithBitSet r = (calcBitSet zeroBits . map fst . VU.toList $ mutationSet r, r)
    -- calculate the bitset pattern for this mutation
    calcBitSet bs [] = bs
    calcBitSet bs (x':xs) =
      let x = maybe (error $ "calcBitSet") id $ B.lookupL mutbit x'
      in  calcBitSet (bs `setBit` x) xs
    -- bijection between mutation position and bit position
    -- @BitSet Bit   <->   Mutated Bit@
    mutbit :: B.BimapHashMap Int Int
    mutbit = B.fromList
           . zipWith (flip (,)) [0 :: Int ..]
           . catMaybes $ zipWith3 genB (BS.unpack origin) (BS.unpack mutation) [0 :: Int ..]
    genB a b k | a == b    = Nothing
               | otherwise = Just $ k

-- | Write a generated landscape to disk.

toFile :: FilePath -> Landscape -> IO ()
toFile fp = BSL.writeFile fp . compress . encodeLazy

fromFile :: FilePath -> IO Landscape
fromFile fp = do
  i <- BSL.readFile fp
  case (decodeLazy . decompress $ i) of
    Left err -> error $ "BioInf.MutationOrder.RNA.fromFile: " ++ err
    Right ls -> return ls

