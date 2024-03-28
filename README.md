# SUL1_natural_variants
Scripts used for analyzing natural variant competition of SUL1. Files in this repository were generated for "High-throughput functional analysis of natural variants in yeast."


mapCCS_reads - aligns PacBio reads to reference sequences. The output files are then ready to be used for barcode-variant mapping. Scripts for doing so can be found here: https://github.com/dunhamlab/PacRAT. To determine how many reads map to reference sequences, first generate reference sequences of the 1,011 collection. generate_fasta_reference

merge_count_barseq - counts reads from barcode sequencing, including trimming and pairing scripts. These outputs were used to calculate fitness.

match_strain_plus_fitness - once barcode-variant map is generated, we can merge fitness and variant data together with strain metadata.
