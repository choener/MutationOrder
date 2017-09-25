
-- | Here we collect the necessary data structures for the RNAs to be
-- compared. This data is serialized to disk once calculated, since it is
-- most likely the part that takes longest.
--
-- TODO if the vienna wrapper allows, we should parallelize the
-- calculations.
--
-- TODO nice interface counting up?

module BioInf.MutationOrder.RNA where

import           Data.Char (toUpper)
import           Codec.Compression.GZip (compress,decompress)
import           Control.Arrow (second)
import           Control.DeepSeq
import           Control.Error
import           Control.Monad (unless)
import           Control.Parallel.Strategies
import           Data.Aeson as DA
import           Data.Bits
import           Data.ByteString (ByteString)
import           Data.Char (isDigit)
import           Data.List (sort,nub,(\\))
import           Data.Maybe (catMaybes)
import           Data.Monoid
import           Data.Serialize
import           Data.Serialize.Instances
import           Data.Text.Encoding (decodeUtf8, encodeUtf8)
import           Data.Tuple (swap)
import           Data.Vector.Serialize
import           Data.Vector.Strategies
import           Debug.Trace
import           GHC.Generics
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy as BSL
import qualified Data.HashMap.Strict as HM
import qualified Data.IntMap.Strict as IM
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import           System.IO.Unsafe (unsafePerformIO)

import qualified Data.Bijection.HashMap as B
import           BioInf.ViennaRNA.Bindings
import qualified Data.PrimitiveArray as PA
import           Biobase.Secondary.Diagrams (D1Secondary(..), mkD1S)



-- | A single RNA with pre-calculated elements.
--
-- All calculations are done at 37 C.
--
-- TODO include the basepair probability matrix? Can we "compress" that
-- one?
--
-- We do not encode D1S into the json

data RNA = RNA
  { mutationSet       :: !(VU.Vector (Int,Char))
    -- ^ we store just the mutation set, since this is more sparse and
    -- gives access to the mutational events.
  , primarySequence   :: !ByteString
    -- ^ store RNA sequence too, for now
  , mfeStructure      :: !ByteString
    -- ^ the mfe structure we get
  , mfeD1S            :: !D1Secondary
    -- ^ efficient structure encoding
  , mfeEnergy         :: !Double
    -- ^ mfe energy of the structure
  , centroidStructure :: !ByteString
    -- ^ the centroid structure
  , centroidD1S       :: !D1Secondary
    -- ^ efficient centroid structure encoding
  , centroidEnergy    :: !Double
  }
  deriving (Show,Eq,Generic)

instance NFData     RNA
instance Serialize  RNA

instance ToJSON RNA where
  toJSON RNA{..} =
    object [ "mutationSet"        .= mutationSet
           , "primarySequence"    .= decodeUtf8 primarySequence
           , "mfeStructure"       .= decodeUtf8 mfeStructure
           , "mfeEnergy"          .= mfeEnergy
           , "centroidStructure"  .= decodeUtf8 centroidStructure
           , "centroidEnergy"     .= centroidEnergy
           ]
  toEncoding RNA{..} =
    pairs (  "mutationSet"        .= mutationSet
          <> "primarySequence"    .= decodeUtf8 primarySequence
          <> "mfeStructure"       .= decodeUtf8 mfeStructure
          <> "mfeEnergy"          .= mfeEnergy
          <> "centroidStructure"  .= decodeUtf8 centroidStructure
          <> "centroidEnergy"     .= centroidEnergy
          )

instance FromJSON RNA where
  parseJSON (Object v) = do
    mutationSet <- v .: "mutationSet"
    primarySequence <- encodeUtf8 <$> v .: "primarySequence"
    let (e,s) = second BS.pack . unsafePerformIO . mfeTemp 37 $ BS.unpack primarySequence
    mfeStructure <- (fmap encodeUtf8 <$> v .:? "mfeStructure") .!= s
    mfeEnergy <- v .:? "mfeEnergy" .!= e
    let (ce,cs) = second BS.pack . unsafePerformIO . centroidTemp 37 $ BS.unpack primarySequence
    centroidStructure <- (fmap encodeUtf8 <$> v .:? "centroidStructure") .!= cs
    centroidEnergy <- v .:? "centroidEnergy" .!= ce
    let mfeD1S = bldD1S mfeStructure
    let centroidD1S = bldD1S centroidStructure
    return RNA{..}

bldD1S :: ByteString -> D1Secondary
bldD1S x = mkD1S (["()"::String], BS.unpack x)

-- | Given the primary sequence and the mutation set, fill the 'RNA'
-- structure.
--
-- NOTE This wraps some @ViennaRNA-bindings@ calls that are in @IO@.
--
-- TODO check if these calls are *really* thread-safe!

mkRNA
  :: Maybe (HM.HashMap ByteString QLine)
  -> ByteString
  -- ^ primary sequence of the *origin* RNA
  -> VU.Vector (Int,Char)
  -- ^ set of mutations compared to the origin
  -> RNA
mkRNA lkup inp' ms = RNA
  { mutationSet       = ms
  , primarySequence   = inp
  , mfeStructure      = mS
  , mfeEnergy         = mE
  , centroidStructure = cS
  , centroidEnergy    = cE
  , mfeD1S            = bldD1S mS
  , centroidD1S       = bldD1S cS
  }
  where
    inp   = insertMutations ms inp'
    ((mE,mS),(cE,cS)) = maybe calculateHere lookup lkup
    calculateHere = ( second BS.pack . unsafePerformIO . mfeTemp 37 $ BS.unpack inp
                    , second BS.pack . unsafePerformIO . centroidTemp 37 $ BS.unpack inp
                    )
    lookup lkup =
      case HM.lookup inp lkup of
        Nothing -> traceShow ("WARNING! have RNA lookup table but have to calculate!", inp) calculateHere
        Just QLine{..} -> (swap qlmfe,swap qlcentroid)

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
  , mutationPositions     :: !(B.BimapHashMap Int Int)
  }
  deriving (Show,Eq,Generic)

instance NFData     Landscape
instance Serialize  Landscape

instance ToJSON Landscape where
  toJSON Landscape{..} =
    object [ "rnas"                 .= rnas
           , "mutationCount"        .= mutationCount
           , "landscapeOrigin"      .= decodeUtf8 landscapeOrigin
           , "landscapeDestination" .= decodeUtf8 landscapeDestination
           , "mutationPositions"    .= mutationPositions
           ]
  toEncoding Landscape{..} =
    pairs (  "rnas"                 .= rnas
          <> "mutationCount"        .= mutationCount
          <> "landscapeOrigin"      .= decodeUtf8 landscapeOrigin
          <> "landscapeDestination" .= decodeUtf8 landscapeDestination
          <> "mutationPositions"    .= mutationPositions
          )

instance FromJSON Landscape where
  parseJSON (Object v) = do
    rnas                  <- v .: "rnas"
    mutationCount         <- v .: "mutationCount"
    landscapeOrigin       <- encodeUtf8 <$> v .: "landscapeOrigin"
    landscapeDestination  <- encodeUtf8 <$> v .: "landscapeDestination"
    mutationPositions     <- v .: "mutationPositions"
    return Landscape{..}

newtype Ancestral = Ancestral ByteString

newtype Extant = Extant ByteString

newtype GlobalBackmutations = GlobalBackmutations Int

newtype BackmutationCol = BackmutationCol Int

-- | The new method for creating the RNA landscape. The outgoing stream of
-- sequences should be directly to the 'BioInf.MutationOrder.SequenceDB'
-- module, and 'writeSequenceFiles' in particular.
--
-- TODO generalise over the alphabet
--
-- TODO currently, the number of backmutations is either @0@ or @1@. For more,
-- a number of things need to be fixed: (i) creation of all variants of
-- @1,2,3...g@ backmutations. That is, consider all sequences that have at
-- least one, but no more than @g@ global back mutations.
--
-- TODO this will create duplicate sequences for those columns where the
-- ancestral and extant have different nucleotides during global backmutation
-- generation.

createRNAlandscape2 ∷ (Monad m) ⇒ [Char] → GlobalBackmutations → [BackmutationCol] → Ancestral → Extant → ExceptT String m (Int,[ByteString])
createRNAlandscape2 alphabet (GlobalBackmutations g) bs (Ancestral a) (Extant e) = do
  -- some sanity checks
  unless (BS.length a == BS.length e) $ throwE "different sequence lengths encountered"
  -- expand later!
  unless (g <= 1) $ throwE "we currently allow *at most* one globally active backmutation"
  -- collect the possible characters for each position.
  let ahm = IM.fromList . zip [0∷Int ..] . map ((:[]) . toUpper) $ BS.unpack a
  let ehm = IM.fromList . zip [0..] . map ( (:[]) . toUpper) $ BS.unpack e
  -- back mutation columns are active together with the above
  let bhm = IM.fromListWith (++) [ (k,alphabet) | BackmutationCol k ← bs ]
  let merge x y = sort . nub $ x++y
  let hm = IM.unionWith merge (IM.unionWith merge ahm ehm) bhm
  -- global back mutations, this will introduce only those characters not already present in each column
  let gbm = IM.fromList [ (k, alphabet \\ hm IM.! k) | k ← [0..BS.length a -1] ]
  -- begin with the set of sequences without any global backmutations
  let localList = map BS.pack . sequence . map snd . IM.toAscList
  let localCount = product . map (length . snd) $ IM.toAscList hm -- do *not* count explicitly!
  let globalModifiers = concat . map (\(k,xs) → map (k,) xs) $ IM.toAscList gbm
  -- we now repeat the local generation, but "splice in" the globally modified characters
  -- TODO we currently allow exactly one global modifier in @[(k,z)]@
  let globals = if g == 1 then [ BS.take k orig `BS.append` (BS.cons z (BS.drop (k+1) orig)) | orig ← localList hm, (k,z) ← globalModifiers ] else []
  let globalCount = if g == 1 then localCount * length globalModifiers else 0
  return $ (localCount + globalCount, localList hm ++ globals)

-- |
--
-- TODO prime candidate for parallelization. ViennaRNA-bindings currently
-- does not allow parallel runs! It would be possible to consider
-- externalizing this, but for now we just run single-threaded.

createRNAlandscape :: Maybe (HM.HashMap ByteString QLine) -> Bool -> ByteString -> ByteString -> (Landscape, [(Int,ByteString)])
createRNAlandscape lkup verbose origin mutation = (ls, zipWith (\mm k -> (k,insertMutations mm origin)) mus [0..])
  where
    ls = Landscape
          { rnas                  = rs -- `using` (parVector chunkSize)
          , mutationCount         = length . filter (>1) . map length $ pms
          , landscapeOrigin       = origin
          , landscapeDestination  = mutation
          , mutationPositions     = mutbit
          }
    rs  = HM.fromList . map pairWithBitSet $ zipWith talk mus [0..]
    talk s c = (if (c `mod` 1000 == 0 && verbose) then traceShow c else id) mkRNA lkup origin s
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
    mutbit = B.fromList
           . zipWith (flip (,)) [0 :: Int ..]
           . catMaybes $ zipWith3 genB (BS.unpack origin) (BS.unpack mutation) [0 :: Int ..]
    genB a b k | a == b    = Nothing
               | otherwise = Just $ k

-- | Write a generated landscape to disk.

toFile :: FilePath -> Landscape -> IO ()
toFile fp = BSL.writeFile fp . compress . encodeLazy

toFileJSON :: FilePath -> Landscape -> IO ()
toFileJSON fp = BSL.writeFile fp . compress . DA.encode

fromFile :: FilePath -> IO Landscape
fromFile fp = (decodeLazy . decompress) <$> BSL.readFile fp >>= \case
  Left err -> error $ "BioInf.MutationOrder.RNA.fromFile: " ++ err
  Right ls -> return ls

fromFileJSON :: FilePath -> IO Landscape
fromFileJSON fp = (DA.eitherDecode' . decompress) <$> BSL.readFile fp >>= \case
  Left err -> error $ "BioInf.MutationOrder.RNA.fromFile: " ++ err
  Right ls -> return ls


-- stupid parsing for quintuple rnafold lines

data QLine = QLine
  { qlSequence  :: ByteString
  , qlmfe       :: (ByteString,Double)
  , qlensemble  :: (ByteString,Double)
  , qlcentroid  :: (ByteString,Double)
  }
  deriving (Show)



qlines f = do
  ls <- BS.lines <$> BS.readFile f
  return $ qlhm $ go ls
  where go [] = []
        go ls = let (hs,ts) = splitAt 5 ls
                in  parseql hs : go ts
        parseql [s,m,e,c,_] =
              QLine s
                    (stupid m)
                    (stupid e)
                    (stupid c)
        stupid bs = let (h:ts) = BS.words bs
                        r = BS.dropWhile (\c -> not $ isDigit c || c=='-') $ BS.unwords ts
                    in (h, read $ BS.unpack $ BS.takeWhile (\c -> isDigit c || c =='-' || c=='.') r)
        qlhm xs = HM.fromList $ map (\q -> (qlSequence q, q)) xs
