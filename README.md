[![Build Status](https://travis-ci.org/choener/MutationOrder.svg?branch=master)](https://travis-ci.org/choener/MutationOrder)

Determine the most likely order of mutations from one RNA sequence to another.

    Walter Costa, Maria Beatriz and Hoener zu Siederdissen, Christian and Tulpan, Dan and Stadler, Peter F. and Nowick, Katja  
    *Uncovering the Structural Evolution of the Human Accelerated Region 1*  
    2017, submitted  
    [preprint](http://www.bioinf.uni-leipzig.de/~choener/pdfs/wal-hoe-2017.pdf)  

# General information

Given two RNA sequences, one ancestral, and one extant, we want to determine
the most likely path of evolution under different measures of fitness.

This program produces the (i) maximum-likelihood path, (ii) all end
probabilities, (iii) all start-end probabilities, (iv) all edge probabilities,
and (v) the maximum expected accuracy path for these two RNA sequences.

In detail:  
(i)   gives the optimal path(s) for the fitness function  
(ii)  gives for each nucleotide polymorphism, how likely it is, that this mutation was introduced *last*  
(iii) looks at all pairs of (first mutation, last mutation) and gives the probability that these two mutations are the begin and end of the chain of mutations  
(iv)  yields for all pairs of nodes (i -> j) the probability that this path occurs, over the whole ensemble of all possible paths  
(v)   produces the path of maximal weight using the probabilities produced in (iv)  

# Usage instructions

## example usage

We assume that you have two Fasta files, *chimp_118.fa* and *human_118.fa* but
they can be named however is convenient. Each file has to contain exactly one
sequence and both sequences have to be of the same length.

For testing with chimp and human, the provided chimp-human.json.gz database
should be used, otherwise the initial foldings will be recalculated. All
required files are available under 'Binaries' at the bottom of the page.

In case, you don't want or can't use the provided work database, run
./MutationOrder with --verbose

We then run

    ./MutationOrder --workdb chimp-human.json.gz --scoretype pairdistcen --onlypositive --outputprefix test chimp_118.fa human_118.fa

This will generate ```test.run```, ```test-edge.eps```, and
```test-meaorder.eps```.

The ```test.run``` file provides extensive output of the optimal path, the
first-last probabilities, the edge probabilities, and the mea output. This
conforms to (i) -- (v) mentioned above.

The two ```eps``` files give a graphical representation of the edge
probabilities, for the ```meaorder``` in order of the path of maximum expected
accuracy.

The work database collects intermediate structures and their folding and is
only created once. The initial run will, however, take some time. I.e. for
'HAR1' this requires 1-4 hours depending on the machine. Further runs complete
*much* faster. In minutes for HAR1.

## Command-line options

    --help        provides short help
    --verbose     will show folding steps during the initial run

    -w
     --workdb=ITEM              the database where to store intermediate foldings
    -t
    --temperature=NUM           annealing temperature. Values close to 0 favor optimal paths. The default is 1.0
    --fillweight=FILLWEIGHT     provides logarithmic and linear fill styles for probability plots. The full style always fills the box
    --fillstyle=FILLSTYLE       normally, boxes are sized, but all in the same color. This changes the opacity of the color as well. Does not work well for eps files
    --cooptcount=INT            how many co-optimals to count for (the count in the .run file is produced differently)
    --cooptprint=INT            how many co-optimals to actually print
    --outprefix=ITEM            how to prefix all output files
    --scoretype=SCORETYPE       choose 'mfe', 'centroid', 'pairdistmfe', or 'pairdistcen' for the evaluation of each mutational step
    --positivesquared           square positive energies to penalize bad moves
    --onlypositive              minimize only over penalties, not energy gains
    --equalstart                each possible mutation is selected with equal probability as the initial one
    --posscaled=NUM,NUM         in =x,y will exponentiate all numbers >=x by the constant y. For value k>=x, we have k^y
    --lkupfile=ITEM             developer option to feed the initial work database with known foldings (usable but very raw and undocumented. needs 5-line rnafold output)
    --showmanual                will show this manual

The allowed score types are:  

    mfe
which optimizes based on the minimum free energy of each intermediate sequence
    centroid
which instead looks at the energy of the centroid structure
    pairdistmfe
which minimizes the base pair distance between following mutations using mfe structures
    pairdistcen
which minimizes the base pair distance between following mutations using centroid structures



# Installation

Follow [this
link](http://www.bioinf.uni-leipzig.de/~choener/software/MutationOrder.html) to
the bottom of the page. Binaries are available for download and installation
from sources via *Haskell Stack* are described.


#### Contact

Christian Hoener zu Siederdissen  
Leipzig University, Leipzig, Germany  
choener@bioinf.uni-leipzig.de  
http://www.bioinf.uni-leipzig.de/~choener/  

