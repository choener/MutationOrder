[![Build Status](https://travis-ci.org/choener/MutationOrder.svg?branch=master)](https://travis-ci.org/choener/MutationOrder)

Determine the most likely order of mutations from one RNA sequence to another.

1.  Walter Costa, Maria Beatriz and Hoener zu Siederdissen, Christian and Tulpan, Dan and Stadler, Peter F. and Nowick, Katja  
    *Uncovering the Structural Evolution of the Human Accelerated Region 1*  
    2017, submitted  
    [preprint](http://www.bioinf.uni-leipzig.de/~choener/pdfs/wal-hoe-2017.pdf)  

# Usage instructions

We assume that you have two Fasta files, *from.fa* and *to.fa* but they can be
named however is convenient (say *chimp.fa* and *human.fa*).  Each file has to
contain exactly one sequence and both sequences have to be of the same length.

We then run

```./MutationOrder --workdb from-to.json.gz --scoretype pairdistcen --onlypositive --outputprefix test```

This will generate ```test.run```, ```test-edge.eps```, and
```test-meaorder.eps```. The ```test.run``` file provides extensive output of
the optimal path, the first-last probabilities, the edge probabilities, and the
mea output. The two ```eps``` files give a graphical representation of the edge
probabilities, for the ```meaorder``` in order of the path of maximum expected
accuracy.

The ```--scoretype``` allows for ```mfe```, ```centroid```, ```pairdistcen```,
and ```pairdistmfe```, which analyse possible evoluationary paths according to
mfe energy, centroid energy, smallest base pair distances for each step in the
```cen```troid or ```mfe``` case.



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

