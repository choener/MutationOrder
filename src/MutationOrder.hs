
{-# Options_GHC -fno-cse #-}

module Main where

import System.Console.CmdArgs
import System.FilePath
import qualified Data.ByteString.Char8 as BS
import Control.Monad

import BioInf.MutationOrder
import BioInf.MutationOrder.RNA (createRNAlandscape)

data ScoreType
  = Mfe
  | Centroid
  | PairDistMfe
  | PairDistCen
  deriving (Show,Data,Typeable)

data Options
  = Options
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
    }
  | GenSequences
    { infiles :: [FilePath]
    }
  deriving (Show,Data,Typeable)

oOptions = Options
  { infiles       = def &= args
  , workdb        = "work.db" &= help "name of the database to store intermediates in"
  , temperature   = 1.0  &= help "lower temperatures favor the more optimal paths, defaults to 1.0"
  , fillweight    = FWlog
  , fillstyle     = FSfull
  , cooptcount    = 100000
  , cooptprint    = 2
  , outprefix     = "tmp"
  , scoretype     = Centroid &= help "choose 'mfe', 'centroid', 'pairdistmfe', or 'pairdistcen' for the evaluation of each mutational step"
  , positivesquared = False &= help "square positive energies to penalize worse structures"
  , onlypositive  = False &= help "minimize only over penalties, not energy gains"
  , equalStart    = False
  , posscaled     = Nothing
  , lkupfile = Nothing
  }

oGenSequences = GenSequences
  { infiles = def &= args
  }

main :: IO ()
main = do
  o <- cmdArgs $ modes [oOptions, oGenSequences] &= verbosity
  case o of
    Options{} -> mainProgram o
    GenSequences{} -> genSequences o

genSequences o = do
  let GenSequences{..} = o
  ancestral <- stupidReader $ infiles !! 0
  current   <- stupidReader $ infiles !! 1
  let ls = snd $ createRNAlandscape Nothing False ancestral current
  forM_ ls $ \(k,sq) -> BS.putStrLn sq
  return ()

mainProgram oOptions = do
  let Options{..} = oOptions
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

