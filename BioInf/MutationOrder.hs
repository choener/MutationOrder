
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

import           Control.Arrow (first,second)
import           Control.Error
import           Control.Monad.IO.Class (liftIO)
import           Control.Monad (unless,forM_,when,forM)
import           Control.Lens
import           Data.ByteString.Strict.Lens
import           Data.Bits
import           Data.ByteString (ByteString)
import           Data.Char (toUpper)
import           Data.Function (on)
import           Data.List (groupBy,sortBy,foldl',(\\),sort)
import           Data.List.Split (chunksOf)
import           Data.Ord (comparing)
import           Data.Tuple (swap)
import           Debug.Trace
import           Numeric.Log
import qualified Control.Parallel.Strategies as Par
import qualified Data.ByteString.Char8 as BS
import qualified Data.HashMap.Strict as HM
import qualified Data.Map.Strict as M
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Data.Trie as Trie
import qualified Data.Vector.Unboxed as VU
import           System.Directory (doesFileExist)
import           System.Exit (exitFailure)
import           System.Exit (exitSuccess)
import           System.IO (withFile,IOMode(WriteMode),hPutStrLn,Handle)
import           Text.Printf

import           ADP.Fusion.Term.Edge.Type (From(..),To(..))
import           Biobase.Secondary.Diagrams (d1Distance)
import           Data.PrimitiveArray (fromEdgeBoundaryFst, EdgeBoundary(..), (:.)(..), getBoundary)
import           Data.PrimitiveArray.ScoreMatrix
import           Diagrams.TwoD.ProbabilityGrid
import qualified Data.Bijection.HashMap as B
import qualified Data.PrimitiveArray as PA
import qualified ShortestPath.SHP.Edge.MinDist as SHP
import           Data.Bits.Ordered

import           BioInf.MutationOrder.EdgeProb
import           BioInf.MutationOrder.MinDist
import           BioInf.MutationOrder.RNA
import           BioInf.MutationOrder.SequenceDB
import qualified BioInf.MutationOrder.BackMutations as BM



runMutationOrder verbose fw fs fwdScaleFunction probScaleFunction cooptCount cooptPrint lkupFile outprefix workdb temperature equalStart [ancestralFP,currentFP] = do
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
    let eplStartWeight = if equalStart
          then Prelude.map (const 1) $ M.elems rowMarginals
          else Prelude.map (Exp . log) $ M.elems rowMarginals
    let eprobsLast = edgeProbScoreMatrix ls eplStartWeight eps
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

-- | This function will run a subset of the possible backmutations. In
-- particular, only those mutations for one particular backmutation column are
-- being used.
--
-- TODO less explicit transformer stack!

runBackmutationVariants
  ∷ BM.ScaleFunction Double
  → FilePath
  -- ^ where the work db lives
  → [Char]
  -- ^
  → Ancestral
  -- ^ The ancestral sequence from which we mutate away
  → Extant
  -- ^ The sequence to which to mutate to
  → Int
  -- ^ The backmutation / intermediate mutation we look at. Indexed with @[1..sequence length]@.
  -- TODO wrap in newtype that enforces this. We have some index structure saying "Start at 1".
  → ExceptT String IO ()
runBackmutationVariants scaleFun workdb alphabet ancestral extant ipos' = do
  let ipos = ipos' - 1
  -- error guarding
  unless (ipos >=0) $ throwE "ipos can not be negative"
  unless (ipos < (BS.length $ getAncestral ancestral)) $ throwE "ipos larger than sequence length"
  unless (BS.length (getAncestral ancestral) == BS.length (getExtant extant)) $ throwE "ancestral and extant sequence do not have equal length"
  -- Load all sequences for the original problem -- they are needed anyway
  -- TODO i think @origSeqs'@ does not hold the original sequences, need to check tonight
  (seqCount, origSeqs', variants) ← createRNAlandscape2 alphabet (Right [ipos]) [] ancestral extant
  let origSeqs = Trie.fromList [ (s,()) | s ← origSeqs' ]
  varSeqs ← forM (alphabet \\ [getAncestral ancestral `BS.index` ipos, getExtant extant `BS.index` ipos]) $ \v → do
    let vsqs = Trie.fromList [ (s,()) | (i,c,ss) ← variants, c == v, s ← ss ]
    return (v, vsqs)
  liftIO $ print workdb
  -- TODO do only a single pass over the data
  -- origStrs ← liftIO $ filter (\r → rnaFoldSequence r `Trie.member` origSeqs) <$> map rna2dna <$> readRNAfoldFiles workdb
  allStrs ← liftIO $ filter (\r → rnaFoldSequence r `Trie.member` origSeqs
                                  || or [rnaFoldSequence r `Trie.member` vs | (_,vs) ← varSeqs] )
          <$> map rna2dna <$> readRNAfoldFiles workdb
  let (rnas,_,_,_) = genSet ancestral extant Nothing allStrs
  -- The @ipos@ declares how many variants we have
  forM_ (alphabet \\ [getAncestral ancestral `BS.index` ipos, getExtant extant `BS.index` ipos]) $ \v → do
    let varSeqs = Trie.fromList [ (s,()) | (i,c,ss) ← variants, c == v, s ← ss ]
    varStrs ← liftIO $ filter (\r → rnaFoldSequence r `Trie.member` varSeqs) <$> map rna2dna <$> readRNAfoldFiles workdb
    let (ntrs,iposbitset,numMuts,mutpos) = genSet ancestral extant (Just (ipos,v)) allStrs
    --liftIO $ print (numMuts, iposbitset, rnas, ntrs, HM.size rnas, HM.size ntrs)
    let fmd = BM.forwardMinDist numMuts
                                scaleFun
                                iposbitset
                                rnas
                                ntrs
    let bts = BM.backtrackMinDist1 numMuts
                                   scaleFun
                                   iposbitset
                                   ipos'
                                   rnas
                                   ntrs
                                   mutpos
                                   fmd
    let evi = BM.forwardEvidence numMuts
                                 (partfun' scaleFun)
                                 iposbitset
                                 rnas
                                 ntrs
    liftIO $ printf "Unobserved Mutation:   Position %3d Nucleotide %c   Delta: %5.1f   lnZ: %10.2f\n"
              ipos' v
              (BM.forwardMinDistValue fmd)
              (ln $ BM.forwardEvidenceValue evi)
    liftIO $ mapM_ T.putStrLn $ take 1 bts
    return ()
  return ()

mfeDelta' :: Bool → Bool → BM.ScaleFunction Double
mfeDelta' mx sq frna trna = (if sq then (\z → if z > 0 then z^2 else z) else id) . (if mx then max 0 else id) $ t - f
  where t = rnaFoldMFEEner trna
        f = rnaFoldMFEEner frna
{-# Inlinable mfeDelta' #-}

centroidDelta' :: Bool → Bool → BM.ScaleFunction Double
centroidDelta' mx sq frna trna = (if sq then (\z → if z > 0 then z^2 else z) else id) . (if mx then max 0 else id) $ t - f
  where t = rnaFoldCentroidEner trna
        f = rnaFoldCentroidEner frna
{-# Inlinable centroidDelta' #-}

partfun' ∷ BM.ScaleFunction Double → BM.ScaleFunction (Log Double)
partfun' f frna trna = Exp . negate $ f frna trna
{-# Inlineable partfun' #-}

rna2dna ∷ RNAfoldResult → RNAfoldResult
rna2dna r = r { rnaFoldSequence = BS.map go $ rnaFoldSequence r } where
  go x =let x' = toUpper x in if x' == 'U' then 'T' else x'

-- | Given ancestral and extant sequence, and possibly an intermediate
-- mutation, as well as a list of intermediates create the set to RNA structure
-- mapping.
--
-- TODO some of the functions here should go into @Lib-SequencePolymorphism@.
--
-- TODO should run within @ExceptT@ !

genSet
  ∷ Ancestral
  → Extant
  → Maybe (Int,Char)
  -- ^ If @Just@ then the 0-based position and character of the intermediate
  -- mutation.
  → [RNAfoldResult]
  → (HM.HashMap Int RNAfoldResult, Int, Int, B.BimapHashMap Int Int)
genSet (Ancestral a') (Extant e') v xs = (HM.fromList kv, ipos, B.size posbit, posbit)
  where kv = [ (b, maybe (kvErr b) id $ HM.lookup (pat2str b) lkupRes) | b ← bits ]
        kvErr b = error $ show (lkupRes, sort $ map rnaFoldSequence xs, v, posbit, b, bits, pat2str b)
        -- update a/e based on if we have the intermediate mutation set up.
        a = maybe a' (\(i,c) → unpackedChars.ix i .~ c $ a') v
        e = maybe e' (\(i,c) → unpackedChars.ix i .~ c $ e') v
        -- all positions where the two bytestrings differ, together with the
        -- differing characters
        ks = filter (\(_,i,j) → i/=j) $ zip3 [0∷Int ..] (BS.unpack a') (BS.unpack e')
        -- turn into a bijection of actual position (first) and bit in bitset
        -- (second)
        posbit ∷ B.BimapHashMap Int Int
        posbit = B.fromList $ zip (ks^..traverse._1) [0∷Int ..]
        -- all bit patterns
        bits = [0 .. 2^B.size posbit - 1]
        -- convert a bit pattern to an actual string, to be looked up. Start
        -- with the ancestral sequence and for each @1@, modify the character
        -- into the one encountered in the extant sequence.
        pat2str ∷ Int → ByteString
        pat2str = let go s k = unpackedChars.ix (lk k) .~ (e `BS.index` lk k) $ s
                      lk = maybe (error "lk") id . B.lookupR posbit
                  in  foldl' go a . activeBitsL
        -- lookup from sequence to RNAfoldResult
        lkupRes = HM.fromList [ (rnaFoldSequence x,x) | x ← xs ]
        -- is the global mutation intermediate? Yes: then >= 0 is the bitset
        -- element to be returned here
        ipos = maybe (-1) (\(z,_) → if BS.index a' z /= BS.index e' z then (maybe (error "ipos") id $ B.lookupL posbit z) else -1) v

-- | Run the intermediate / backmutation order variant. This variant is slow,
-- and requires large pre-calculated files, we parallelize and aggregate as
-- much as possible.
--
-- TODO read monad ?!

{-
runBackmutationVariants
  ∷ Int
  → [Char]
  → GlobalBackmutations
  → [BackmutationCol]
  → Ancestral
  → Extant
  → ExceptT String IO ()
runBackmutationVariants aggregate alphabet globback backcols ancestral extant = do
  -- Load all sequences for the original problem -- they are needed anyway
  (seqCount, origSeqs', variants) ← createRNAlandscape2 alphabet globback backcols ancestral extant
  let origSeqs = Trie.fromList [ (s,()) | s ← origSeqs' ]
  origStrs ← filter (\r → rnaFoldSequence r `Trie.member` origSeqs) <$> readRNAfoldFiles (error "workdb")
  let rnas = error "bitset -> rnafoldresult data"
  -- Group into sets of @aggregate@ elements for sequence aggregation
  let ass = chunksOf aggregate variants
  forM_ ass $ \as → do
    -- the required sequences are given by @origSeqs@ but modified at the appropriate position
    let aggrSeqss = map (\(p,n,xs) → (p,n,Trie.fromList [ (x,()) | x ← xs ])) as
    let allss = foldl' (\z (_,_,x) → Trie.unionL z x) Trie.empty aggrSeqss
    -- read in the structures for all as
    aggrStr ← filter (\r → rnaFoldSequence r `Trie.member` allss) <$> readRNAfoldFiles (error "workdb")
    let go (p,n,xs) = (p,n) where
          -- prepare the @ntrs@ data structure for each as
          ntrs = undefined $ filter (\r → rnaFoldSequence r `Trie.member` xs) aggrStr
          fMnD = BM.forwardMinDist (error "number of known mutations") (error "scale function") rnas ntrs
          fwdZ = BM.forwardEvidence (error "number of known mutations") (error "scale function for evidence") rnas ntrs
    -- parallel calculation and output for the different cases
    let rs = Par.parMap Par.rdeepseq go aggrSeqss
    forM_ rs $ \r → do
      return ()
    return ()
  -- print each output
  return ()
-}

