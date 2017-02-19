
-- | Run all steps of the HoxCluster algorithms in order.
--
-- This will produce the following:
-- 
-- 1. run the minimal distance algorithm, give the minimal distance score
-- and return all co-optimal paths
--
-- 2. run the end-probability algorithm and return the probability that
-- each node is the begin/end of a chain
--
-- 3. run the edge probability algorithm and give the probability for each
-- @from :-> to@ edge
--
-- 4. with the edge probabilities, run the maximal probability path
-- algorithm, return that probability and all co-optimal paths
--
-- TODO -Pretty should yield a structure to be given to the eps or svg
-- generator. This allows more flexibility. Does diagrams offer
-- serialization?
--
-- TODO All this should be wrapped and available as a function. not just
-- providing output files.

module BioInf.MutationOrder
  ( module BioInf.MutationOrder
  , FillWeight (..)
  , FillStyle (..)
  , ScaleFunction (..)
  ) where

import           Numeric.Log
import           Data.List (groupBy,sortBy)
import           Data.Function (on)
import           Control.Monad (unless,forM_)
import           Data.ByteString (ByteString)
import qualified Data.ByteString.Char8 as BS
import           System.Directory (doesFileExist)
import           System.Exit (exitFailure)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import           Text.Printf
import           Data.Ord (comparing)
import qualified Data.Map.Strict as M

import qualified Data.Bijection.HashMap as B
import           Diagrams.TwoD.ProbabilityGrid
import           Data.PrimitiveArray (fromEdgeBoundaryFst, EdgeBoundary(..))
import qualified ShortestPath.SHP.Edge.MinDist as SHP

import           BioInf.MutationOrder.RNA
import           BioInf.MutationOrder.MinDist
import           BioInf.MutationOrder.EdgeProb



runMutationOrder verbose fw fs scaleFunction cooptCount cooptPrint workdb temperature [ancestralFP,currentFP] = do
  --
  -- Initial stuff and debug information
  --
  ancestral <- stupidReader ancestralFP
  current   <- stupidReader currentFP
  ls <- withDumpFile workdb ancestral current $ createRNAlandscape verbose ancestral current
  print $ mutationCount ls
  --
  -- Run co-optimal lowest energy changes
  --
  let (e,bs) = runCoOptDist scaleFunction ls
  printf "Best energy gain: %10.4f\n" e
  printf "Number of co-optimal paths: %10d\n" (length $ take cooptCount bs)
  putStrLn ""
  forM_ (take cooptPrint bs) T.putStrLn
  --
  -- Run edge probability Inside/Outside calculations. These take quite
  -- a while longer.
  --
  let (ibs,eps) = edgeProbPartFun scaleFunction temperature ls
  let mpks = sortBy (comparing snd) . B.toList $ mutationPositions ls
  putStr "       "
  forM_ mpks $ \(mp,k) -> printf " %6d" k
  putStrLn ""
  putStr "       "
  forM_ mpks $ \(mp,k) -> printf " %6d" (mp+1)
  putStrLn ""
  forM_ (zip (groupBy ((==) `on` (fromEdgeBoundaryFst . fst)) eps) mpks) $ \(rps,(mp,k)) -> do
    let (eb,_) = head rps
    printf "%3d %3d" k (mp+1)
    forM_ rps $ \(eb,Exp p) -> printf (" %6.4f") (exp p)
    printf "   %6.4f" (Prelude.sum $ map (exp . ln . snd) rps)
    printf "\n"
  let colSums = M.fromListWith (+) [ (c,p) | ((_ :-> c),p) <- eps ]
  putStr "    Σ  "
  forM_ (M.toList colSums) $ \(c,Exp p) -> printf (" %6.4f") (exp p)
  putStrLn "\n"
  --
  -- Generate the path with maximal edge probability
  --
  let eprobs = edgeProbScoreMatrix ls eps
  let (Exp maxprob,mpbt) = SHP.runMaxEdgeProb eprobs
  printf "Maximal Edge Log-Probability Sum: %6.4f\n" maxprob
  putStrLn "From extant to ancestral, most recent to first mutation:\n"
  forM_ mpbt $ \bt -> T.putStrLn bt
  putStrLn ""
{-# NoInline runMutationOrder #-}

-- | Stupid fasta reader

stupidReader :: FilePath -> IO ByteString
stupidReader fp = do
  inp <- BS.lines <$> BS.readFile fp
  let xs = filter (\x -> not (BS.null x) && BS.head x /= '>') inp
  return $ BS.concat xs

-- | @withDumpFile@ is like @idIO :: a -> IO a@ in that it returns the data
-- we give to the function. However, in case the dump file exists, we read
-- it and return its contents, instead of recalculating. If it does not
-- exist, we dump the data in addition to returning it. This forces the
-- @Landscape@.

withDumpFile
  :: FilePath
  -- ^ The path we store the serialized and compressed dump in
  -> ByteString
  -- ^ ancestral / origin sequence
  -> ByteString
  -- ^ destination sequence
  -> Landscape
  -- ^ the element which is to be serialized in the dump, or which would be
  -- the data in the dump
  -> IO Landscape
  -- ^ the data we put in, but maybe taken from the dump file
withDumpFile fp ancestral current l = do
  dfe <- doesFileExist fp
  if dfe then do
    printf "using database %s to load sequence information\n" fp
    ls <- fromFileJSON fp
    -- now we check if we have a sane DB file
    unless (landscapeOrigin ls == ancestral && landscapeDestination ls == current) $ do
      putStrLn "ancestral or target sequence do not match those stored in the work database"
      putStrLn $ "given ancestral: " ++ BS.unpack ancestral
      putStrLn $ "DB    ancestral: " ++ (BS.unpack $ landscapeOrigin ls)
      putStrLn $ "given current:   " ++ BS.unpack current
      putStrLn $ "DB    current:   " ++ (BS.unpack $ landscapeDestination ls)
      exitFailure
    return ls
  else do
    printf "database %s does not exist! Folding all intermediate structures. This may take a while!\n" fp
    toFileJSON fp l
    return l

