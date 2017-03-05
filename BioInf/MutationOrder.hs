
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

import qualified Data.Vector.Unboxed as VU
import           Data.Tuple (swap)
import           Control.Monad (unless,forM_,when)
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
import           Control.Arrow (first,second)
import           System.IO (withFile,IOMode(WriteMode),hPutStrLn,Handle)
import           System.Exit (exitSuccess)

import           ADP.Fusion.Term.Edge.Type (From(..),To(..))
import           Data.PrimitiveArray (fromEdgeBoundaryFst, EdgeBoundary(..), (:.)(..), getBoundary)
import           Data.PrimitiveArray.ScoreMatrix
import qualified Data.PrimitiveArray as PA
import           Diagrams.TwoD.ProbabilityGrid
import qualified Data.Bijection.HashMap as B
import qualified ShortestPath.SHP.Edge.MinDist as SHP
import           Biobase.Secondary.Diagrams (d1Distance)

import           BioInf.MutationOrder.EdgeProb
import           BioInf.MutationOrder.MinDist
import           BioInf.MutationOrder.RNA



runMutationOrder verbose fw fs fwdScaleFunction probScaleFunction cooptCount cooptPrint lkupFile outprefix workdb temperature [ancestralFP,currentFP] = do
  -- only run if out file(s) do not exist
  dfe <- doesFileExist (outprefix ++ ".run")
  when dfe $ do
    printf "%s.run exists, ending now!\n" outprefix
    exitSuccess
  withFile (outprefix ++ ".run") WriteMode $ \oH -> do
    printf "%s.run job started!\n" outprefix
    --
    -- Initial stuff and debug information
    --
    ancestral <- stupidReader ancestralFP
    current   <- stupidReader currentFP
    lkup <- case lkupFile of {Nothing -> return Nothing; Just f -> Just <$> qlines f}
    ls <- withDumpFile oH workdb ancestral current . fst $ createRNAlandscape lkup verbose ancestral current
    let mpks = sortBy (comparing snd) . B.toList $ mutationPositions ls
    let bitToNuc = M.fromList $ map (swap . first (+1)) mpks
    let nn = length mpks
    hPrintf oH "number of mutations: %d\n" $ mutationCount ls
    hPrintf oH "\n%s\n\n" $ replicate 80 '='
    --
    -- Run co-optimal lowest energy changes
    --
    let (e,bs) = runCoOptDist fwdScaleFunction ls
    let (ecount,countcount) = runCount fwdScaleFunction ls
    -- split co-optimals into "want to print" and "want to count";
    -- @countbs@ should be possible to stream
    let (printbs,countbs) = splitAt cooptPrint bs
    hPrintf oH "Best energy gain: %10.4f\n" e
    hPrintf oH "Number of co-optimal paths: %10d\n" countcount -- ((length printbs) + (length $ take (cooptCount-cooptPrint) bs))
    forM_ printbs (T.hPutStrLn oH)
    hPrintf oH "%s\n\n" $ replicate 80 '='
    --
    -- Run @First@ probability algorithm to determine the probability for
    -- each mutation to be the initial one
    --
    -- TODO this is completely wrong, because it still starts at the
    -- ancestral sequence. We would have to start at the extant sequence.
    -- Need to later think about this. But do not use any @First@ functions
    -- now!
    {-
    hPrintf oH "Chain begin probabilities:\n"
    let fps = boundaryPartFunFirst Nothing probScaleFunction ls
    forM_ mpks $ \(mp,k) -> hPrintf oH "  %6d" (mp+1)
    hPrintf oH "\n"
    forM_ fps $ \(_, Exp p) -> hPrintf oH "  %6.4f" (exp p)
    hPrintf oH "\n\n"
    printf "\n"
    -}
    --
    -- Run @Last@ probability algorithm to determine the probability for
    -- each mutation to be the last one
    --
    hPrintf oH "Chain end probabilities:\n"
    let fps = boundaryPartFunLast Nothing probScaleFunction ls
    forM_ mpks $ \(mp,k) -> hPrintf oH "  %6d" (mp+1)
    hPrintf oH "\n"
    forM_ (bpNormalized fps) $ \(_, Exp p) -> hPrintf oH "  %6.4f" (exp p)
    hPrintf oH "\n\n%s\n\n" $ replicate 80 '='
    --printf "\n"
    --
    -- Run specialized versions of the above, restricting the first mutation
    -- to the given one. Marginalized over the last probability, and rescaled
    -- we get the first probability. Completely printed out, we get the joint
    -- probability for each @i,j@ to be @first,last@ in the chain.
    --
    hPrintf oH "Restricted chain end probabilities\n"
    let rbps = map (\(mp,k) -> (mp,k,boundaryPartFunLast (Just k) probScaleFunction ls)) mpks
    {-
    forM_ rbps $ \(mp,k,bp) -> do
      hPrintf oH "%5d %5d\n" (mp+1) k
      forM_ (bpUnnormalized bp) $ \(l,Exp p) -> hPrintf oH "%7d " (bitToNuc M.! getBoundary l)
      hPrintf oH "\n"
      forM_ (bpUnnormalized bp) $ \(l,p) -> hPrintf oH "%7.2f " (exp . ln $ p / bpTotal bp)
      hPrintf oH "\n"
    hPrintf oH "\n"
    -}
    -- collect all restricted partition function scores and prepare for
    -- normalization
    let firstlastUn = M.fromList [ ((mp+1,bitToNuc M.! getBoundary l), logp)
                                 | (mp,k,bp) <- rbps, (l,logp) <- bpUnnormalized bp
                                 ]
    let firstlastZ = Numeric.Log.sum [ bpTotal bp | (_,_,bp) <- rbps ]
    let firstlastLogP = M.map (/firstlastZ) firstlastUn
    let firstlastP = M.map (exp . ln) firstlastLogP
    -- rowMarginals gives the total probability that the mutation order
    -- begins with this mutation.
    let rowMarginals = M.mapKeysWith (+) fst firstlastP
    -- colMarginals gives the total probability that the mutation order
    -- ends with this mutation.
    let colMarginals = M.mapKeysWith (+) snd firstlastP
    hPrintf oH "       "
    forM_ (M.elems bitToNuc) $ \mut -> hPrintf oH "%6d " mut
    hPrintf oH "         Σ\n"
    forM_ (M.elems bitToNuc) $ \frst -> do
      hPrintf oH "%4d   " frst
      forM_ (M.elems bitToNuc) $ \lst -> hPrintf oH "%6.4f " (firstlastP M.! (frst,lst))
      hPrintf oH "    %6.4f\n" $ rowMarginals M.! frst
    hPrintf oH "Σ      "
    forM_ (M.elems colMarginals) $ hPrintf oH "%6.4f "
    hPrintf oH "\n\n"
    hPrintf oH "divergence from proper normalization: %10.8f\n" (1 - Prelude.sum firstlastP)
    hPrintf oH "row marginal sum %10.8f\n" (Prelude.sum rowMarginals)
    hPrintf oH "col marginal sum %10.8f\n" (Prelude.sum colMarginals)
    hPrintf oH "\n%s\n\n" $ replicate 80 '='
    -- debug on
    {-
    hPrintf oH "%f\n" $ ln firstlastZ
    hPrintf oH "%s " $ replicate 10 ' '
    forM_ (M.elems bitToNuc) $ \mut -> hPrintf oH "%10d " mut
    hPrintf oH "\n"
    forM_ (M.elems bitToNuc) $ \frst -> do
      hPrintf oH "%8d   " frst
      forM_ (M.elems bitToNuc) $ \lst -> hPrintf oH "%10.4f " (ln $ firstlastUn M.! (frst,lst))
      hPrintf oH "\n"
    hPrintf oH "\n"
    hPrintf oH "%f\n" $ ln firstlastZ
    hPrintf oH "%s " $ replicate 10 ' '
    forM_ (M.elems bitToNuc) $ \mut -> hPrintf oH "%10d " mut
    hPrintf oH "\n"
    forM_ (M.elems bitToNuc) $ \frst -> do
      hPrintf oH "%8d   " frst
      forM_ (M.elems bitToNuc) $ \lst -> hPrintf oH "%10.4f " ((ln $ firstlastUn M.! (frst,lst)) - ln firstlastZ)
      hPrintf oH "\n"
    hPrintf oH "\n"
    -}
    -- debug off
    -- debug on
    -- calculate first weight, unnormalized
  --  let firstUn = M.fromList [ ]
    -- debug off
    --
    --
    -- Run edge probability Inside/Outside calculations. These take quite
    -- a while longer.
    --
    let (ibs,eps) = edgeProbPartFun probScaleFunction ls
    hPrintf oH "pairwise next mutation probabilities:\n\n"
    hPrintf oH "       "
    forM_ mpks $ \(mp,k) -> hPrintf oH " %6d" k
    hPutStrLn oH ""
    hPrintf oH "       "
    forM_ mpks $ \(mp,k) -> hPrintf oH " %6d" (mp+1)
    hPutStrLn oH ""
    forM_ (zip (groupBy ((==) `on` (fromEdgeBoundaryFst . fst)) eps) mpks) $ \(rps,(mp,k)) -> do
      let (eb,_) = head rps
      hPrintf oH "%3d %3d" k (mp+1)
      forM_ rps $ \(eb,Exp p) -> hPrintf oH (" %6.4f") (exp p)
      hPrintf oH "   %6.4f" (Prelude.sum $ map (exp . ln . snd) rps)
      hPrintf oH "\n"
    let colSums = M.fromListWith (+) [ (c,p) | ((_ :-> c),p) <- eps ]
    hPrintf oH "    Σ  "
    forM_ (M.toList colSums) $ \(c,Exp p) -> hPrintf oH (" %6.4f") (exp p)
    hPutStrLn oH "\n"
    gridFile [SVG,EPS] (outprefix ++ "-edge") fw fs nn nn (map (show . (+1) . fst) mpks) (map (show . (+1) . fst) mpks) (map snd eps)
    hPrintf oH "\n%s\n\n" $ replicate 80 '='
    --
    -- Generate the path with maximal edge probability
    --
    {-
    let eprobsFirst = edgeProbScoreMatrix ls (Prelude.map (Exp . log) $ M.elems colMarginals) eps
    let (Exp maxprob,mpbt) = SHP.runMaxEdgeProbFirst eprobsFirst
    hPrintf oH "Maximal Edge Log-Probability Sum: %6.4f with at least %d co-optimal paths\n" maxprob (length $ take cooptCount mpbt)
    hPutStrLn oH "first mutation to extant species\n"
    forM_ (take cooptPrint mpbt) $ \bt -> do
      let extractMut (SHP.BTnode (_:.To n)) = n
          extractMut (SHP.BTedge (From ff:.To tt)) = ff
      let mutationOrder = tail $ scanl (\set mut -> set `setBit` extractMut mut) zeroBits (reverse bt)
      let prettyPrint mut k = do
            let rna = rnas ls HM.! mut
            hPrintf oH "   %3s  %s\n        %s   MFE %6.4f\n        %s   CNT %6.4f\n"
                    (maybe "anc" (show . (+1) . fst . (!!) mpks) k)
                    (BS.unpack $ primarySequence rna)
                    (BS.unpack $ mfeStructure rna)
                    (mfeEnergy rna)
                    (BS.unpack $ centroidStructure rna)
                    (centroidEnergy rna)
            hPutStrLn oH $ replicate 8 ' ' ++ (take (BS.length $ primarySequence rna) . concat $ zipWith (\xs x -> xs ++ show x) (repeat $ "    .    ") (drop 1 $ cycle [0..9]))
      prettyPrint zeroBits Nothing
      forM_ (zip (reverse bt) mutationOrder) $ \case
        (SHP.BTnode (_:.To n),mut) -> prettyPrint mut (Just n)
        (SHP.BTedge (From ff:.To tt),mut) -> prettyPrint mut (Just ff)
      hPutStrLn oH ""
    hPutStrLn oH ""
    -}
    -- the rowMarginals hold the probabily to begin with a mutation. Since
    -- @Last@ goes from first to last mutation, this is what we need.
    let eprobsLast = edgeProbScoreMatrix ls (Prelude.map (Exp . log) $ M.elems rowMarginals) eps
    --print eprobsLast
    --print $ PA.assocs $ scoreMatrix eprobsLast
    let (Exp maxprobLast,lastLogProbs,mpbtLast') = SHP.runMaxEdgeProbLast eprobsLast
    let mpbtLast = map reverse mpbtLast'
    --print maxprobLast
    --print lastLogProbs
    --mapM_ print $ mpbtLast
    hPrintf oH "Fraction of optimal choice for each final mutation:\n"
    forM_ lastLogProbs $ \(PA.Boundary b, _) -> hPrintf oH "  %6d" $ bitToNuc M.! b
    hPrintf oH "\n"
    forM_ lastLogProbs $ \(_, p) -> hPrintf oH "  %6.4f" $ exp $ ln (p / Exp maxprobLast)
    hPrintf oH "\n\n"
    hPrintf oH "Maximal edge log-probability sum: %6.4f (P = %10.8f) with at least %d co-optimal paths\n" maxprobLast (exp maxprobLast) (length $ take cooptCount mpbtLast)
    hPutStrLn oH "(first mutation to extant species)\n"
    forM_ (take cooptPrint mpbtLast) $ \bt -> do
      let extractMut (SHP.BTnode (_:.To n)) = n
          extractMut (SHP.BTedge (From ff:.To tt)) = tt
      let mutationOrder = tail $ scanl (\set mut -> set `setBit` extractMut mut) zeroBits bt
      let prettyPrint mut k = do
            let rna = rnas ls HM.! mut
            hPrintf oH "   %3s  %s\n        %s   MFE %6.4f\n        %s   CNT %6.4f\n"
                    (maybe "anc" (show . (+1) . fst . (!!) mpks) k)
                    (BS.unpack $ primarySequence rna)
                    (BS.unpack $ mfeStructure rna)
                    (mfeEnergy rna)
                    (BS.unpack $ centroidStructure rna)
                    (centroidEnergy rna)
            hPutStrLn oH $ replicate 8 ' ' ++ (take (BS.length $ primarySequence rna) . concat $ zipWith (\xs x -> xs ++ show x) (repeat $ "    .    ") (drop 1 $ cycle [0..9]))
      prettyPrint zeroBits Nothing
      forM_ (zip bt mutationOrder) $ \case
        (SHP.BTnode (_:.To n),mut) -> prettyPrint mut (Just n)
        (SHP.BTedge (From ff:.To tt),mut) -> prettyPrint mut (Just tt)
      hPutStrLn oH ""
    hPutStrLn oH ""
    let meaOrder =
          let go = \case SHP.BTnode (_:.To n) -> n
                         SHP.BTedge (From ff:.To tt) -> tt
          in  map go $ concat $ take 1 mpbtLast
    let meaAnno = map (\k -> map (show . (+1) . fst) mpks !! k) meaOrder
    let meaEps = [ (ee !! k) !! l | let ee = groupBy ((==) `on` (fromEdgeBoundaryFst . fst)) eps, k <- meaOrder, l <- meaOrder ]
    gridFile [SVG,EPS] (outprefix ++ "-edge-meaorder") fw fs nn nn meaAnno meaAnno (map snd meaEps)
    --print eps
    --print meaEps
    {-
    let eprobsLast = edgeProbScoreMatrix ls (Prelude.map (Exp . log) $ M.elems rowMarginals) eps
    let (Exp maxprobLast,lastLogProbs,mpbtLast) = SHP.runMaxEdgeProbLast eprobsLast
    print maxprobLast
    print $ map (\(k,Exp p) -> (k,exp $ p - maxprobLast)) lastLogProbs
    mapM_ print $ concat $ take 2 mpbtLast
    let eprobsLast = edgeProbScoreMatrix ls (Prelude.map (Exp . log) $ M.elems colMarginals) eps
    let (Exp maxprobLast,lastLogProbs,mpbtLast) = SHP.runMaxEdgeProbLast eprobsLast
    print maxprobLast
    print $ map (\(k,Exp p) -> (k,exp $ p - maxprobLast)) lastLogProbs
    mapM_ print $ concat $ take 2 mpbtLast
    -}
{-# NoInline runMutationOrder #-}

posScaled :: Double -> Double -> ScaleFunction -> ScaleFunction
posScaled l s = scaleByFunction go where
  go d | d >= l    = d ** s
       | otherwise = d
  {-# Inline go #-}
{-# Inlinable posScaled #-}

-- | Basepair distance

basepairDistanceMFE :: ScaleFunction
basepairDistanceMFE frna trna = fromIntegral $ d1Distance (mfeD1S frna) (mfeD1S trna)

basepairDistanceCentroid :: ScaleFunction
basepairDistanceCentroid frna trna = fromIntegral $ d1Distance (centroidD1S frna) (centroidD1S trna)

-- | Scale function for normal mfe delta energies

mfeDelta :: ScaleFunction
mfeDelta frna trna = mfeEnergy trna - mfeEnergy frna
{-# Inlinable mfeDelta #-}

-- | Scale function for normal centroid delta energies

centroidDelta :: ScaleFunction
centroidDelta frna trna = centroidEnergy trna - centroidEnergy frna
{-# Inlinable centroidDelta #-}

-- | Square positive "contributions", making bad moves more unlikely

squaredPositive :: ScaleFunction -> ScaleFunction
squaredPositive sf = scaleByFunction sp sf where
  sp d
    | d > 0     = d * d
    | otherwise = d
  {-# Inline sp #-}
{-# Inlinable squaredPositive #-}

-- | Scale by temperature (for probability stuff)

scaleTemperature :: Double -> ScaleFunction -> ScaleFunction
scaleTemperature t sf = scaleByFunction (/t) sf
{-# Inlinable scaleTemperature #-}

scaleByFunction f sf = \frna trna ->
  let d = sf frna trna
  in  f d
{-# Inlinable scaleByFunction #-}

-- | Basepair distance

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
  :: Handle
  -> FilePath
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
withDumpFile oH fp ancestral current l = do
  dfe <- doesFileExist fp
  if dfe then do
    hPrintf oH "using database %s to load sequence information\n" fp
    ls <- fromFileJSON fp
    -- now we check if we have a sane DB file
    unless (landscapeOrigin ls == ancestral && landscapeDestination ls == current) $ do
      hPutStrLn oH "ancestral or target sequence do not match those stored in the work database"
      hPutStrLn oH $ "given ancestral: " ++ BS.unpack ancestral
      hPutStrLn oH $ "DB    ancestral: " ++ (BS.unpack $ landscapeOrigin ls)
      hPutStrLn oH $ "given current:   " ++ BS.unpack current
      hPutStrLn oH $ "DB    current:   " ++ (BS.unpack $ landscapeDestination ls)
      exitFailure
    return ls
  else do
    hPrintf oH "database %s does not exist! Folding all intermediate structures. This may take a while!\n" fp
    toFileJSON fp l
    return l

