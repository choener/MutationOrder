
module Main where

import System.Console.CmdArgs
import System.FilePath

import BioInf.MutationOrder

data ScoreType
  = Mfe
  | Centroid
  | StructureDistance
  deriving (Show,Data,Typeable)

data Options = Options
  { infiles       :: [FilePath]
  , workdb        :: FilePath
  , temperature   :: Double
  , fillweight    :: FillWeight
  , fillstyle     :: FillStyle
  , cooptcount    :: Int
  , cooptprint    :: Int
  , figurenames   :: String
  , scoretype     :: ScoreType
  , positivesquared :: Bool
  , onlypositive  :: Bool
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
  , figurenames   = "fig-"
  , scoretype     = Centroid &= help "choose 'mfe', 'centroid', or 'structuredistance' for the evaluation of each mutational step"
  , positivesquared = False &= help "square positive energies to penalize worse structures"
  , onlypositive  = False &= help "minimize only over penalties, not energy gains"
  } &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs oOptions
  isL <- isLoud
  let fwdScaleFunction
        = (if positivesquared then squaredPositive else id)
        . (if onlypositive then (scaleByFunction (max 0)) else id)
        $ (case scoretype of Mfe -> mfeDelta; Centroid -> centroidDelta;)
  let insideScaleFunction
        = scaleTemperature temperature
        . (if positivesquared then squaredPositive else id)
        . (if onlypositive then (scaleByFunction (max 0)) else id)
        $ (case scoretype of Mfe -> mfeDelta; Centroid -> centroidDelta;)
  runMutationOrder isL fillweight fillstyle fwdScaleFunction insideScaleFunction cooptcount cooptprint figurenames workdb temperature infiles

