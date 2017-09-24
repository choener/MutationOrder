
-- | Functions to deal with a (large) RNA sequence data base and corresponding
-- RNAfold foldings.

module BioInf.MutationOrder.SequenceDB where

import           Codec.Compression.GZip (compress)
import           Control.Error
import           Control.Monad.IO.Class (liftIO)
import           Data.ByteString.Char8 (ByteString)
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as BSL
import           System.Directory (doesFileExist)
import           System.FilePath ((</>), (<.>))



type PrefixLen = Int
type SeqsPerFile = Int

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

