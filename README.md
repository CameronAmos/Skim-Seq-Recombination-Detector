# Skim-Seq-Recombination-Detector

A script used to identify recombination breakpoints using parental variant markers of skim-seq sequence data derived from a single inbred line.
The script identifies gross transition points within the data to identify both homozygous and heterozygous regions accurately.
Transitions in data are identified and scored, with the sharpest transitions scoring the highest. 
A logic tree then filters these candidate recombination sites, filtering out lower-scoring duplicates and impossible recombinations based on surrounding data.
If a recombiantion cannot be found between two high-confident regions, a second search ensues that finds the most favorable transition point
using a large sliding window, its size based on the proximity of more confident recombination breakpoints.
It is highly configurable, and works with high accuracy on skim-seq data only 0.028x in coverage.

Input file example containing the chromosome, position, and if it was the Reference or Alternative allele.
chr1A_TA299   323292  P2
chr1A_TA299   323424  P2
chr1A_TA299   323528  P2
chr1A_TA299   323867  P1
...    ... ...

Output file 1 example of the xy data, containing the chromosome, marker index, and roaming score (useful for plotting the function)

chr1A_TA299     0       1

chr1A_TA299     1       2

chr1A_TA299     2       3

chr1A_TA299     3       4

...    ... ...

Output file 2 example of the detected recombinations containing:}

Chromosome, left marker pos, right marker pos, average pos, transition type, marker index, and breakpoint score (out of 2)

chr1A_TA299     529372358       529428788       529400573       P1/P2   109847  1.963

chr2A_TA299     54816056        54820632        54818344        P2/P1   4710    1.927

chr3A_TA299     42509727        42527879        42518803        P1/P2   5213    1.968

chr4A_TA299     60813562        60824007        60818784        P2/P1   6995    1.939

...    ... ...
