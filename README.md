# multiple-barcoded-seqs-demultiplexing
Assign CRISPR data into bins, given an amplicons set.

In order to demultiplex each of the libraries we have applied the following workflow:

1. Align library reads against all 136 sequences of interest (88 MRE amplicons + 48 calibration reads) using the Needleman-Wunsch algorithm. The Needleman-Wunsch algorithm was the ideal option since, as a global alignment technique, it can capture alignments having potentially big gaps (which we expect to have due to CRISPR/Cas9 activity).

2. Following alignment, each read is either assigned to a specific MRE bin or is classified as calibration read.

3. The barcode of each calibration read is checked (4 or 3 nt long) and the read is assigned to a specific calibration bin.

4. For each MRE/calibration, reads with a number of mismatches over 10 nt are discarded since they correspond to unsuccessful alignments.

5. All filtered reads of each bin are aligned with BLAST against the respective wild type sequence of the bin in order to discard any remaining sequencing artefacts.


Step 5 of the demultiplexing process assures that mis-classified reads in step 1 (that are most likely sequencing artefacts) are actually discarded (Fig. 1.2) and the clean depth is retrieved for each MRE, normalised by read depth across all libraries (Fig 1.3).
