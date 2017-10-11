
{-# Options_GHC -fno-cse #-}

module Main where

import           Control.Monad.IO.Class (liftIO, MonadIO)
import           Control.Monad
import           Control.Error
import           Data.FileEmbed
import qualified Data.ByteString.Char8 as BS
import           System.Console.CmdArgs
import           System.Directory (createDirectoryIfMissing)
import           System.Exit (exitSuccess, exitFailure)
import           System.FilePath

import BioInf.MutationOrder
import BioInf.MutationOrder.RNA
import BioInf.MutationOrder.SequenceDB

data Options
  = MutationOrder
    { infiles       :: [FilePath]
    , workdb        :: FilePath
    , temperature   :: Double
    , fillweight    :: FillWeight
    , fillstyle     :: FillStyle
    , cooptcount    :: Int
    , cooptprint    :: Int
    , outprefix     :: FilePath
    , scoretype     :: ScoreType
    , positivesquared :: Bool
    , onlypositive  :: Bool
    , equalStart    :: Bool
    , posscaled :: Maybe (Double,Double)
    , lkupfile :: Maybe FilePath
    , showmanual    :: Bool
    }
  | GenSequences
    { ancestralSequence ∷ FilePath
    -- ^ the (presumably) ancestral sequence from which to start the mutation
    -- order prediction.
    , extantSequence ∷ FilePath
    -- ^ the extant or target sequence to which to mutate.
    -- 
    , globalbackmutations ∷ Int
    -- ^ how many global back mutations to allow. For example, if set to @1@,
    -- then for each position in the sequence all "unobserved" nucleotides are
    -- possible with one switch to and one switch back. It is generally a bad
    -- idea to have more than @1@ here and we only allow exactly 1 backmutation
    -- with a specialized algorithm (that still takes quite a while to run).
    , backmutationcolumns ∷ [Int]
    -- ^ Additional columns for backmutations, may overlap with observed
    -- mutations.
    , sequenceLimit ∷ Int
    -- ^ Complain if the number of sequences is above this limit.
    , workdb ∷ FilePath
    -- ^ will write files to @workdb </> seqs@
    , prefixlength ∷ Int
    -- ^ how many characters as prefix for writing
    , seqsperfile ∷ Int
    -- ^ how many sequences to write into each file
    , alphabet ∷ String
    -- ^ the different characters that are allowed
    }
  | Backmutation
    { ancestralSequence ∷ FilePath
    -- ^ the (presumably) ancestral sequence from which to start the mutation
    -- order prediction.
    , extantSequence ∷ FilePath
    -- ^ the extant or target sequence to which to mutate.
    , workdb ∷ FilePath
    -- ^
    , position ∷ Int
    -- ^
    , alphabet ∷ String
    , scoretype     :: ScoreType
    , positivesquared :: Bool
    , onlypositive  :: Bool
    }
  deriving (Show,Data,Typeable)

oMutationOrder = MutationOrder
  { infiles       = def &= args
  , workdb        = "work.db" &= help "name of the database to store intermediates in"
  , temperature   = 1.0  &= help "lower temperatures favor the more optimal paths, defaults to 1.0"
  , fillweight    = FWlog &= help "scale options for probability plots: fwlog, fwlinear, fwfill"
  , fillstyle     = FSfull &= help "fill options for probability plots: fsopacitylog, fsopacitylinear, fsfull"
  , cooptcount    = 1000  &= help "how many co-optimals to count"
  , cooptprint    = 2   &= help "how many co-optimals to actually print out"
  , outprefix     = "tmp" &= help "prefix for output files"
  , scoretype     = Centroid &= help "choose 'mfe', 'centroid', 'pairdistmfe', or 'pairdistcen' for the evaluation of each mutational step"
  , positivesquared = False &= help "square positive energies to penalize worse structures"
  , onlypositive  = False &= help "minimize only over penalties, not energy gains"
  , equalStart    = False &= help "run mea with equal start probabilities"
  , posscaled     = Nothing &= help "--posscaled=x,y   scale all values >= x by using y as exponent"
  , lkupfile = Nothing  &= help "developer option: if an RNAfold file with foldings exists, then use it"
  , showmanual = False  &= help "shows the manual"
  }

oGenSequences = GenSequences
  { ancestralSequence = def
  , extantSequence = def
  , globalbackmutations = 0
  , backmutationcolumns = []
  , sequenceLimit = 1000000
  , workdb = def
  , prefixlength = 4
  , seqsperfile = 10000
  , alphabet = ""
  }

oBackmutation = Backmutation
  { ancestralSequence = def
  , extantSequence = def
  , workdb = def
  , position = -1
  , alphabet = ""
  , scoretype     = Centroid &= help "choose 'mfe', 'centroid', 'pairdistmfe', or 'pairdistcen' for the evaluation of each mutational step"
  , positivesquared = False &= help "square positive energies to penalize worse structures"
  , onlypositive  = False &= help "minimize only over penalties, not energy gains"
  }

main :: IO ()
main = do
  o <- cmdArgs $ modes [oMutationOrder &= auto, oGenSequences, oBackmutation] &= verbosity
  case o of
    MutationOrder{} → mainProgram o
    GenSequences{}  → genSequences o
    Backmutation{}  → runBackmutation o

-- | This is a simple wrapper around the RNA landscape creation. Landscape
-- creation generates all sequences between ancestral and extant sequence. It
-- will take into account additional global backmutations and further
-- backmutation columns. Note that this can *very easily* lead to combinatorial
-- explosion.
--
-- Needs extra options on: (i) globally active backmutations. For each position
-- try all four nucleotides. (ii) Additional active columns. For each position,
-- try all four nucleotides.

genSequences o = do
  let GenSequences{..} = o
  e ← runExceptT $ do
    when (null ancestralSequence) $ throwE "ancestral sequence file?"
    when (null extantSequence) $ throwE "extant sequence file?"
    when (null workdb) $ throwE "work db directory?"
    when (null alphabet) $ throwE "use --alphabet=ACGT (or ACGU if your fasta files are RNA-based)"
    a ← liftIO $ stupidReader ancestralSequence
    e ← liftIO $ stupidReader extantSequence
    (numSeqs, origs, sqs) ← createRNAlandscape2
                              alphabet
                              (Left $ GlobalBackmutations globalbackmutations)
                              (map BackmutationCol backmutationcolumns)
                              (Ancestral a) (Extant e)
    unless (numSeqs <= sequenceLimit) $ throwE $ "combinatiorial explosion (" ++ show numSeqs ++ "): reduce search space or allow for higher --sequencelimit"
    -- write out sequences
    let sdb = workdb </> "sequences"
    liftIO $ createDirectoryIfMissing True sdb
    writeSequenceFiles sdb prefixlength seqsperfile
      $ origs ++ concatMap (\(_,_,xs) → xs) sqs
  case e of
    Left err → print err >> exitFailure
    Right () → return ()
  {-
  ancestral <- stupidReader $ infiles !! 0
  current   <- stupidReader $ infiles !! 1
  let ls = snd $ createRNAlandscape Nothing False ancestral current
  forM_ ls $ \(k,sq) -> BS.putStrLn sq
  return ()
  -}

embeddedManual = $(embedFile "README.md")

mainProgram oOptions = do
  let MutationOrder{..} = oOptions
  when showmanual $ do
    BS.putStrLn embeddedManual
    exitSuccess
  unless (length infiles == 2) $ do
    BS.putStrLn embeddedManual
    putStrLn "\n\n\nThis program expects exactly two equal-length fasta files as input"
    exitFailure
  isL <- isLoud
  runMutationOrder isL fillweight fillstyle scoretype positivesquared posscaled onlypositive cooptcount cooptprint lkupfile outprefix workdb temperature equalStart infiles

runBackmutation ∷ Options → IO ()
runBackmutation Backmutation{..} = do
  e ← runExceptT $ do
    when (null ancestralSequence) $ throwE "ancestral sequence file?"
    when (null extantSequence) $ throwE "extant sequence file?"
    when (null workdb) $ throwE "work db directory?"
    when (null alphabet) $ throwE "use --alphabet=ACGT (or ACGU if your fasta files are RNA-based)"
    a ← liftIO $ Ancestral <$> stupidReader ancestralSequence
    e ← liftIO $ Extant    <$> stupidReader extantSequence
    let scaleFun = case scoretype of
                      Centroid → centroidDelta' onlypositive positivesquared
                      Mfe → mfeDelta' onlypositive positivesquared
                      PairDistMfe → mfebpdist' onlypositive positivesquared
                      PairDistCen → centroidbpdist' onlypositive positivesquared
    runBackmutationVariants scaleFun workdb alphabet a e position
  case e of
    Left err → print (err ∷ String) >> exitFailure
    Right () → return ()
  return ()

