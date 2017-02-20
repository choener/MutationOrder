
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

import           Control.Monad (unless,forM_)
import           Data.Bits
import           Data.ByteString (ByteString)
import           Data.Function (on)
import           Data.List (groupBy,sortBy)
import           Data.Ord (comparing)
import           Numeric.Log
import qualified Data.ByteString.Char8 as BS
import qualified Data.HashMap.Strict as HM
import qualified Data.Map.Strict as M
import qualified Data.Text as T
import qualified Data.Text.IO as T
import           System.Directory (doesFileExist)
import           System.Exit (exitFailure)
import           Text.Printf

import           ADP.Fusion.Term.Edge.Type (From(..),To(..))
import           Data.PrimitiveArray (fromEdgeBoundaryFst, EdgeBoundary(..), (:.)(..))
import           Diagrams.TwoD.ProbabilityGrid
import qualified Data.Bijection.HashMap as B
import qualified ShortestPath.SHP.Edge.MinDist as SHP

import           BioInf.MutationOrder.EdgeProb
import           BioInf.MutationOrder.MinDist
import           BioInf.MutationOrder.RNA



runMutationOrder verbose fw fs scaleFunction cooptCount cooptPrint fignames workdb temperature [ancestralFP,currentFP] = do
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
  let nn = length mpks
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
  gridFile [SVG,EPS] (fignames ++ "edge") fw fs nn nn (map (show . (+1) . fst) mpks) (map (show . (+1) . fst) mpks) (map snd eps)
  --
  -- Generate the path with maximal edge probability
  --
  let eprobs = edgeProbScoreMatrix ls eps
  let (Exp maxprob,mpbt) = SHP.runMaxEdgeProb eprobs
  printf "Maximal Edge Log-Probability Sum: %6.4f with at least %d co-optimal paths\n" maxprob (length $ take cooptCount mpbt)
  putStrLn "first mutation to extant species\n"
  forM_ (take cooptPrint mpbt) $ \bt -> do
    let extractMut (SHP.BTnode (_:.To n)) = n
        extractMut (SHP.BTedge (From ff:.To tt)) = ff
    let mutationOrder = tail $ scanl (\set mut -> set `setBit` extractMut mut) zeroBits (reverse bt)
    let prettyPrint mut k = do
          let rna = rnas ls HM.! mut
          printf "   %3s  %s\n        %s   MFE %6.4f\n        %s   CNT %6.4f\n"
                  (maybe "anc" (show . (+1) . fst . (!!) mpks) k)
                  (BS.unpack $ primarySequence rna)
                  (BS.unpack $ mfeStructure rna)
                  (mfeEnergy rna)
                  (BS.unpack $ centroidStructure rna)
                  (centroidEnergy rna)
    prettyPrint zeroBits Nothing
    forM_ (zip (reverse bt) mutationOrder) $ \case
      (SHP.BTnode (_:.To n),mut) -> prettyPrint mut (Just n)
      (SHP.BTedge (From ff:.To tt),mut) -> prettyPrint mut (Just ff)
    putStrLn ""
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

