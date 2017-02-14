
module Main where

import System.Console.CmdArgs
import System.FilePath

import BioInf.MutationOrder



data Options = Options
  { infiles     :: [FilePath]
  , workdb      :: FilePath
  , temperature :: Double
  , fillweight  :: FillWeight
  , fillstyle   :: FillStyle
  }
  deriving (Show,Data,Typeable)

oOptions = Options
  { infiles     = def &= args
  , workdb      = def &= help "name of the database to store intermediates in"
  , temperature = 0.01  &= help "lower temperatures favor the more optimal paths, defaults to 0.01"
  , fillweight  = FWlog
  , fillstyle   = FSfull
  } &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs oOptions
  isL <- isLoud
  return ()
  runMutationOrder isL fillweight fillstyle workdb temperature infiles

