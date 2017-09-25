
-- | Functions to deal with a (large) RNA sequence data base and corresponding
-- RNAfold foldings.

module BioInf.MutationOrder.SequenceDB where

import           Codec.Compression.GZip (compress,decompress)
import           Control.Error
import           Control.Monad.IO.Class (liftIO)
import           Data.ByteString.Char8 (ByteString)
import           Data.Char (isDigit)
import qualified Data.Char as C
import qualified Data.Attoparsec.ByteString.Lazy as A
import qualified Data.Attoparsec.ByteString.Char8 as AC
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as BSL
import           System.Directory (doesFileExist, getDirectoryContents)
import           System.FilePath ((</>), (<.>))
import           System.IO.Unsafe (unsafeInterleaveIO)



type PrefixLen = Int
type SeqsPerFile = Int

-- * Plain sequences, with no folding yet

-- | Writes sequence files out to disk. Does *not* check if duplicates are
-- written. The user should do this externally via @sort@ and @uniq@. The next
-- free prefix is chosen.

writeSequenceFiles
  ∷ FilePath
  -- ^ Directory to write to.
  → PrefixLen
  -- ^ Number of characters in the prefix to count up. Fail if the prefix space
  -- is too small.
  → SeqsPerFile
  -- ^ How many sequences to write per file
  → [ByteString]
  -- ^ Sequences to write out.
  → ExceptT String IO ()
writeSequenceFiles fp pfx spf xs = do
  -- cull prefix list to unused prefixes only
  let unused [] = throwE "writeSequenceFiles: prefix space empty"
      unused (x:xs) = do
        dfe ← liftIO . doesFileExist $ fp </> x <.> "gz"
        if dfe
          then unused xs
          else return $ x:xs
  -- will write out file names to unused prefixes.
  let go ps [] = return ()
      go ps ys = do
        let (here,there) = splitAt spf ys
        (p:qs) ← unused ps
        -- write to compressed file @p@
        let fname = fp </> p <.> "gz"
        liftIO . BSL.writeFile fname . compress . BSL.unlines $ map BSL.fromStrict here
        go qs there
  let pfxs = sequence $ replicate pfx ['a'..'z']
  go pfxs xs

-- | Reads all sequences.
--
-- TODO make sure to read lazily enough that collection happens without first
-- reading in the whole list.

readSequenceFiles
  ∷ FilePath
  → ExceptT () IO [ByteString]
readSequenceFiles fp = do
  return undefined

-- * Structural information from RNAfold

-- | Space-efficient RNAfold structure

data RNAfoldResult = RNAfoldResult
  { rnaFoldSequence       ∷ !ByteString
  , rnaFoldMFEStruc       ∷ !ByteString
  , rnaFoldMFEEner        ∷ !Double
  , rnaFoldMfeFrequency   ∷ !Double
  , rnaFoldEnsembleStruc  ∷ !ByteString
  , rnaFoldEnsembleEner   ∷ !Double
  , rnaFoldCentroidStruc  ∷ !ByteString
  , rnaFoldCentroidEner   ∷ !Double
  , rnaFoldDiversity      ∷ !Double
  }
  deriving (Read,Show,Eq,Ord)

-- | Lazily read @RNAfold@ structures.
--
-- TODO use @pipes/machines@!
-- TODO use filemanip or something like that

readRNAfoldFiles
  ∷ FilePath
  → ExceptT String IO [RNAfoldResult]
readRNAfoldFiles fp = do
  let go [] = return []
      go (f:fs) = do
        bs ← liftIO $ unsafeInterleaveIO (decompress <$> BSL.readFile f)
        rs ← bslToRNAfoldResult bs
        rss ← go fs
        return $ rs ++ rss
  liftIO (getDirectoryContents fp) >>= go
{-# NoInline readRNAfoldFiles #-}

-- |

bslToRNAfoldResult ∷ (Monad m) ⇒ BSL.ByteString → ExceptT String m [RNAfoldResult]
bslToRNAfoldResult bs = do
  case A.eitherResult $ A.parse pRNAfold bs of
    Left e  → throwE e
    Right r → return r
{-# Inline bslToRNAfoldResult #-}

-- |
--
-- @
-- echo "CCCAAAGGG\nCCCAAAGGG" | ./RNAfold -p
-- CCCAAAGGG
-- (((...))) ( -1.20)
-- (((...))) [ -1.41]
-- (((...))) { -1.20 d=1.06}
--  frequency of mfe structure in ensemble 0.707288; ensemble diversity 1.67  
-- CCCAAAGGG
-- (((...))) ( -1.20)
-- (((...))) [ -1.41]
-- (((...))) { -1.20 d=1.06}
--  frequency of mfe structure in ensemble 0.707288; ensemble diversity 1.67  
-- @

pRNAfold ∷ A.Parser [RNAfoldResult]
pRNAfold = A.many' go where
  go = do
    -- 1. sequence
    rnaFoldSequence       ← AC.takeWhile AC.isAlpha_ascii <* AC.skipSpace
    -- 2. mfe
    rnaFoldMFEStruc       ← AC.takeTill AC.isSpace <* AC.skipSpace
    rnaFoldMFEEner        ← AC.char '(' *> AC.skipSpace *> AC.double <* AC.char ')' <* AC.skipSpace
    -- 3. ensemble
    rnaFoldEnsembleStruc  ← AC.takeTill AC.isSpace <* AC.skipSpace
    rnaFoldEnsembleEner   ← AC.char '[' *> AC.skipSpace *> AC.double <* AC.char ']' <* AC.skipSpace
    -- 4. centroid
    rnaFoldCentroidStruc  ← AC.takeTill AC.isSpace <* AC.skipSpace
    rnaFoldCentroidEner   ← AC.char '{' *> AC.skipSpace *> AC.double <* AC.skipSpace
    dequal                ← AC.string "d=" *> AC.double <* AC.char '}' <* AC.skipSpace
    -- 5.mfe frequency and diversity
    AC.string " frequency of mfe structure in ensemble"
    rnaFoldMfeFrequency ← AC.double
    AC.string "; ensemble diversity"
    rnaFoldDiversity ← AC.double
    AC.skipSpace
    return RNAfoldResult{..}
{-# Inline pRNAfold #-}

