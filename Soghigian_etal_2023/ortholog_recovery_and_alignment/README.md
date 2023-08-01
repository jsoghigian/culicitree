Scripts related to ortholog recovery, alignment, and quality control are stored here.

oma_to_orthograph.txt - A description of loops and scripts necessary to adjust OMA results to Orthograph format.

align_and_trim.bash - A script to align and trim putative orthologs with MAFFT and Trimal.  Also requires the python script fasta_subset.py.

summarize_aa.bash - Create concatenated alignments from the output of align_and_trim.bash. Requires fasta_to_catnex and fastaaa_to_catnex, as well as part_converter.sh.

outlier_detection.R - The R commands used to identify outlier sequences based on our branch-length method.
