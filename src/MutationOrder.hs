
module Main where

import System.Console.CmdArgs
import System.FilePath

import BioInf.MutationOrder



data Options = Options
  { infiles       :: [FilePath]
  , workdb        :: FilePath
  , temperature   :: Double
  , fillweight    :: FillWeight
  , fillstyle     :: FillStyle
  , scalefunction :: ScaleFunction
  , cooptcount    :: Int
  , cooptprint    :: Int
  , figurenames   :: String
  }
  deriving (Show,Data,Typeable)

oOptions = Options
  { infiles       = def &= args
  , workdb        = "work.db" &= help "name of the database to store intermediates in"
  , temperature   = 1.0  &= help "lower temperatures favor the more optimal paths, defaults to 0.01"
  , fillweight    = FWlog
  , fillstyle     = FSfull
  , scalefunction = ScaleId
  , cooptcount    = 100000
  , cooptprint    = 2
  , figurenames   = "fig-"
  } &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs oOptions
  isL <- isLoud
  runMutationOrder isL fillweight fillstyle scalefunction cooptcount cooptprint figurenames workdb temperature infiles

