
{-# Options_GHC -fno-cse #-}

module Main where

import           Control.Monad
import           Data.FileEmbed
import qualified Data.ByteString.Char8 as BS
import           System.Console.CmdArgs
import           System.Exit (exitSuccess, exitFailure)
import           System.FilePath

import BioInf.MutationOrder
import BioInf.MutationOrder.RNA (createRNAlandscape)

data ScoreType
  = Mfe
  | Centroid
  | PairDistMfe
  | PairDistCen
  deriving (Show,Data,Typeable)

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
    { infiles :: [FilePath]
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
  { infiles = def &= args
  }

main :: IO ()
main = do
  o <- cmdArgs $ modes [oMutationOrder &= auto, oGenSequences] &= verbosity
  case o of
    MutationOrder{} -> mainProgram o
    GenSequences{} -> genSequences o

genSequences o = do
  let GenSequences{..} = o
  ancestral <- stupidReader $ infiles !! 0
  current   <- stupidReader $ infiles !! 1
  let ls = snd $ createRNAlandscape Nothing False ancestral current
  forM_ ls $ \(k,sq) -> BS.putStrLn sq
  return ()

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
  let fwdScaleFunction
        = (if positivesquared then squaredPositive else id)
        . (maybe id (uncurry posScaled) posscaled)
        . (if onlypositive then (scaleByFunction (max 0)) else id)
        $ (case scoretype of Mfe -> mfeDelta
                             Centroid -> centroidDelta
                             PairDistMfe -> basepairDistanceMFE
                             PairDistCen -> basepairDistanceCentroid)
  let insideScaleFunction
        = scaleTemperature temperature
        . (if positivesquared then squaredPositive else id)
        . (maybe id (uncurry posScaled) posscaled)
        . (if onlypositive then (scaleByFunction (max 0)) else id)
        $ (case scoretype of Mfe -> mfeDelta
                             Centroid -> centroidDelta
                             PairDistMfe -> basepairDistanceMFE
                             PairDistCen -> basepairDistanceCentroid)
  runMutationOrder isL fillweight fillstyle fwdScaleFunction insideScaleFunction cooptcount cooptprint lkupfile outprefix workdb temperature equalStart infiles

