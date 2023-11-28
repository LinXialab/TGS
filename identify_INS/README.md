## System requirements

`python 3.7.4`

python package `pysam` and `pandas`

## Installation guide

Install pysam via `conda` by

`conda install pysam pandas`

## Example data

Example data is provided to test the somatic INS calling. 

The expected output is a candidate somatic INS txt file.

## Instructions for use

### Candidate somatic INS identification ###

Input is expected in BAM format

### manual check ###

You need to manually reviewed the somatic INSs using the Integrative Genomics Viewer based on the 
split alignments of INSs.

### Local assembly of somatic INS ###
1. the INS-supported reads for each somatic INS were extracted based on the record of Sniffles and were assembled using Shasta (version 0.5.1) with the following parameters:

   ```
    • Reads.minReadLength (1000)
    • MinHash.minHashIterationCount (60)
    • MinHash.minFrequency (2)
    • Align.minAlignedMarkerCount (100)
    • MarkerGraph.minCoverage (5)
    ```
    
2. To further improve the quality of the assembled Shasta sequences, Racon (version 1.4.3)  and Medaka (version 1.0.3) were used to polish the draft assemblies.
