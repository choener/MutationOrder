#!/bin/bash
# this is all --noLP, but allowing closing GU pairs

# we have two data sources, second one is not used anymore
DB="chimp-human"
# chimp-human-noclosinggu"

# we have a number of score types
ST="mfe centroid pairdistmfe pairdistcen"

# run with only positive contributions
OP="false true"

# different temperatures to run on
T="1.0 0.9 0.8 0.5 0.1"

# positive values scaled
PS="1,2"

parallel ./MutationOrder options --outprefix HAR1/runs/{1}-{2}-op{3}-T{4}-psno --workdb HAR1/workdb/{1}.json.gz HAR1/fasta/chimp_118.fa HAR1/fasta/human_118.fa --scoretype={2} --onlypositive={3} --temperature={4} ::: $DB ::: $ST ::: $OP ::: $T

parallel ./MutationOrder options --outprefix HAR1/runs/{1}-{2}-op{3}-T{4}-ps{5} --workdb HAR1/workdb/{1}.json.gz HAR1/fasta/chimp_118.fa HAR1/fasta/human_118.fa --scoretype={2} --onlypositive={3} --temperature={4} --posscaled={5} ::: $DB ::: $ST ::: $OP ::: $T ::: $PS

#./MutationOrder options --workdb HAR1/workdb/chimp-human.json.gz HAR1/fasta/chimp_118.fa HAR1/fasta/human_118.fa --outprefix HAR1/runs/chimp-human-mfe -s mfe
